#include "EventSelector.h"

int plotLLratio_wire(const char* fname, const char* templates_sig, const char* templates_bkg, const char* outfname) {

    EventSelector sel(fname, templates_sig, templates_bkg, outfname);

    sel.Loop();

    return 0;
}

int EventSelector::Init()
{
    //==================================================
    // Output file and Hitograms
    //==================================================
    gStyle->SetPalette(kBlueRedYellow);
    gStyle->SetNumberContours(255);
    gStyle->SetHistFillStyle(0);
    gStyle->SetOptStat(0);

    foutf = new TFile(foutfname,"recreate");

    TH1* hKlike = new TH1F("hKlike", "Log-likelihood ratio background/signal K-like hypthesis (Forward)",
			   200, -100, 100);
    TH1* hPlike = new TH1F("hPlike", "Log-likelihood ratio background/signal P-like hypthesis (Backward)",
			   200, -100, 100);
    TH1* hDLLratio = new TH1F("hDLLratio", "Double Log-likelihood ratio background/signal Forward + Backward",
			   200, -100, 100);

    fHists["hKlike"] = hKlike;
    fHists["hPlike"] = hPlike;
    fHists["hDLLratio"] = hDLLratio;

    fHists["hNewTemplate"] = new TH2D("hNewTemplate",
				      ";Number of Wires from Vertex;dE/dx",
				      50,0,50,50,0,50);
    fHists["hNewTemplate_rev"] = new TH2D("hNewTemplate_rev",
					  ";Number of Wires from Oposite End;dE/dx",
					  50,0,50,50,0,50);

    const char* dir_label[] = { "Forward", "Backward" };
    const char* signal_label[] = { "Signal", "Background" };

    for (int i = 0; i < 2; i++) { // fw vs bw
	for (int j = 0; j < 2; j++) { // signal vs bkg
	    TString histname(Form("hLL%d%d", i, j));
	    cout << histname << endl;
	    fHists[histname] = new TH1I(histname,
					Form("%s %s;Log-Likelihood",
					     dir_label[i], signal_label[j]),
					125, -150., 100.);
	}
    }



    return 0;
}

int EventSelector::ProcessEvent(NDKAna* evt)
{
    // retrieved output histograms
    TH1* hKlike    = fHists["hKlike"];
    TH1* hPlike    = fHists["hPlike"];
    TH1* hDLLratio = fHists["hDLLratio"];

    //==================================================
    // Calculate log-likelihood from templates
    //==================================================
    double like_p = 0;
    double like_k = 0;
    double loglike = 0;
    double like_p_rev = 0;
    double like_k_rev = 0;
    double loglike_rev =0;
    double PID_loglike = 0;
    int npoints = evt->n_cal_points[fShortIdx];
    double track_length = evt->track_length[fShortIdx];

    // unnecessary, calo points always sorted such that end of track comes first
    // int points_from_vertex = (evt->track_range[fShortIdx][0] >
    // 			       evt->track_range[fShortIdx][npoints-1]);


    // compare reconstructed dE/dX vs wire number distance from muon vertex with the templates
    for( int k = 0; k < npoints; k++ ){
	// loop over calorimetric points of the best plane

	if( k > 49 ) break; // templates do only 50 wires max

	// indexing from the muon vertex
	int idx = 0;
	if ( fShortDir ==  -1 )
	    idx = npoints - k - 1 ;
	else
	    idx = k;
	idx = (fCaloOrder == 1) ? idx : npoints - 1 - idx ;


	double dEdx =evt->track_dE_dx[fShortIdx][idx];
	fHists["hNewTemplate"] -> Fill(k, dEdx);

	int bin_p = fHProbP ->GetYaxis()->FindBin(dEdx);
	int bin_k = fHProbK ->GetYaxis()->FindBin(dEdx);

	double tmp_like_p = fHProbP->GetBinContent(k+1, bin_p);
	double tmp_like_k = fHProbK->GetBinContent(k+1, bin_k);
	// deal with calo pts where dedx is out of our range - assign
	// them default minimal likelihood
	if (bin_p == 0 || bin_p > fHProbP ->GetYaxis()->GetNbins() ) {
	    tmp_like_p = -5;
	    tmp_like_k = -5;
	}

	like_p += tmp_like_p;
	like_k += tmp_like_k;

	//using reverse template
	int rev_idx = npoints - idx - 1;
	double dEdx_rev =evt->track_dE_dx[fShortIdx][rev_idx];
	fHists["hNewTemplate_rev"] -> Fill(k, dEdx_rev);
	int bin_p_rev = fHProbP_rev ->GetYaxis()->FindBin(dEdx_rev);
	int bin_k_rev = fHProbK_rev ->GetYaxis()->FindBin(dEdx_rev);

	double tmp_like_p_rev = fHProbP_rev->GetBinContent(k+1, bin_p_rev);
	double tmp_like_k_rev = fHProbK_rev->GetBinContent(k+1, bin_k_rev);
	if (bin_p_rev == 0 || bin_p_rev > fHProbP_rev ->GetYaxis()->GetNbins() ) {
	    tmp_like_p_rev = -5;
	    tmp_like_k_rev = -5;
	}

	like_p_rev += tmp_like_p_rev;
	like_k_rev += tmp_like_k_rev;

	// print out info for events where loglikelihood > 0
	// if (tmp_like_p > 0 || tmp_like_k > 0 || tmp_like_p_rev > 0 || tmp_like_k_rev > 0) {
	//     cout<<setw(15)<<"p-like"<<setw(15)<<tmp_like_p<<endl
	// 	<<setw(15)<<"k-like"<<setw(15)<<tmp_like_k<<endl
	// 	<<setw(15)<<"p-like rev"<<setw(15)<<tmp_like_p_rev<<endl
	// 	<<setw(15)<<"k-like rev"<<setw(15)<<tmp_like_k_rev<<endl;
	//     cout<<setw(15)<<"k"<<setw(15)<<k<<endl
	// 	<<setw(15)<<"idx"<<setw(15)<<idx<<endl
	// 	<<setw(15)<<"rev idx"<<setw(15)<<rev_idx<<endl
	// 	<<setw(15)<<"dEdX"<<setw(15)<<dEdx<<endl
	// 	<<setw(15)<<"rev dEdX"<<setw(15)<<dEdx_rev<<endl;
	//     cout<<"--"<<endl;
	// }

    }
    loglike = like_p-like_k;
    loglike_rev = like_p_rev-like_k_rev;
    PID_loglike = loglike + loglike_rev;

    fHists["hLL00"] -> Fill(like_k);
    fHists["hLL01"] -> Fill(like_p);
    fHists["hLL10"] -> Fill(like_k_rev);
    fHists["hLL11"] -> Fill(like_p_rev);

    hKlike->Fill(loglike);
    hPlike->Fill(loglike_rev);
    hDLLratio->Fill(PID_loglike);

    return 0;
}


int EventSelector::Finalise()
{
    // retrieved output histograms
    TH1* hKlike    = fHists["hKlike"];
    TH1* hPlike    = fHists["hPlike"];
    TH1* hDLLratio = fHists["hDLLratio"];

    hKlike->SetLineColor(kRed);
    hPlike->SetLineColor(kBlue);

    TString name = foutfname;
    name.ReplaceAll(".root", "");

    TCanvas* c = new TCanvas("c","", 600, 300);
    c->Divide(2,1);
    c->cd(1);
    THStack* hs = new THStack("hs", "Log-likelihood ratio Background/Signal");
    hs->Add(hPlike);
    hs->Add(hKlike);
    hs->Draw("nostack");
    TLegend* leg = new TLegend(0.7,0.7,0.85,0.85);
    leg->AddEntry(hKlike, "Forward", "l");
    leg->AddEntry(hPlike, "Backward", "l");
    leg->Draw();

    c->cd(2);
    hDLLratio->Draw();

    c->SaveAs(name+"_likelihood_ratios.pdf");

    c->cd(1);
    fHists["hNewTemplate"]->Draw("colz");
    c->cd(2);
    fHists["hNewTemplate_rev"]->Draw("colz");
    c->SaveAs(name+"_dedx_vs_wire.pdf");

    TCanvas* c2 = new TCanvas("c2","", 600, 600);
    c2->Divide(2,2);
    c2->cd(1);
    fHProbK->Draw("colz");
    c2->cd(2);
    fHProbP->Draw("colz");
    c2->cd(3);
    fHProbK_rev->Draw("colz");
    c2->cd(4);
    fHProbP_rev->Draw("colz");
    c2->SaveAs(name+"_prob_templates.pdf");


    auto hs1 = new THStack("hs1","individual log-likelihoods");
    for (int i = 0; i < 2; i++) { // fw vs bw
	for (int j = 0; j < 2; j++) { // signal vs bkg
	    c2->cd(i + j*2 + 1);
	    TString histname = Form("hLL%d%d", i, j);
	    fHists[histname] -> Draw();
	    hs1->Add(fHists[histname]);
	    fHists[histname]->SetLineColor(i + j*2 + 1);
	}
    }
    c2->SaveAs(name+"_likelihoods.pdf");


    // all in one canvas
    auto c3 = new TCanvas("c3", "",600,900);
    c3 -> Divide (2,3);
    c3 -> cd(1);
    hs1->Draw("nostack");
    leg = gPad->BuildLegend();

    c3 -> cd (3);
    hs -> Draw("nostack");
    leg = gPad->BuildLegend();

    c3 -> cd (4);
    hDLLratio->Draw();

    c3->cd(5);
    fHists["hNewTemplate"]->Draw("colz");
    c3->cd(6);
    fHists["hNewTemplate_rev"]->Draw("colz");

    c3->SaveAs(name+"_all_plots.pdf");



    // write out loglikelihood histograms (from the templates)
    fHProbP->Write();
    fHProbK->Write();
    fHProbP_rev->Write();
    fHProbK_rev->Write();

    foutf->Write();
    foutf->Close();


    return 0;
}
