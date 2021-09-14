#ifndef EVENTDISPLAY_H
#define EVENTDISPLAY_H

#include "NDKAna.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

class EventDisplay
{
public:
    EventDisplay(NDKAna* anaSel);

    int next(const int entry = -1);
    int plot();
    void PrintInfo();
    int GetEntry(int entry);

    void AllowBranches(const std::vector<string>& br_list);

    std::vector<std::vector<TArrow*> > GetShowers();
    TMultiGraph** plotTracks ();
    TMultiGraph** plotTracksWV ();
    TMultiGraph** plotMCparts();
    TMultiGraph** plotMCpartsWV();
    TH1* plot2Dgraphs(vector<TGraph2D*>& grv);

private:
    void reverseYaxix(TMultiGraph* mg);

private:
    NDKAna* mAnaSel;
    TTree* mTree;
    int mEntry;

    TCanvas* mCanvas; // 3D->2D projections
    TCanvas* mCanvasWv; // wire view
    TCanvas* mCanvas3d; // 3D view

    // U and V axis direction
    TVector3 u_dir;
    TVector3 v_dir;
};

EventDisplay::EventDisplay(NDKAna* anaSel) :
    mAnaSel(anaSel),
    mTree(anaSel->fChain),
    mEntry(0),
    mCanvas(0),
    mCanvas3d(0)
{
    (u_dir = TVector3(0,0,1)).RotateX(-35.7/180. * TMath::Pi());
    (v_dir = TVector3(0,0,1)).RotateX(35.7/180. * TMath::Pi());

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    mCanvas = new TCanvas("c", "2D View", 1024,256);
    mCanvas->Divide(4,1);

    mCanvasWv = new TCanvas("cwv", "2D Wire-Time View", 450,500);
    mCanvasWv->Divide(1,3);

    mCanvas3d = new TCanvas("c3d", "3D View", 300,300);

    mTree->SetBranchStatus("*", 0);
    // allow only necessary branches
    vector<string> allowed ({"subRunNo"
		, "eventNo"
		, "n_reco_tracks"
		//, "n_cal_points_byplane"
		//, "track_calo_xyz_byplane"
		, "n_track_points"
		, "track_point_xyz"
		, "track_vtx"
		, "track_end"
		, "mc_npart"
		, "mc_pdg"
		, "mc_id"
		, "mc_npoints"
		, "mc_xyz"
		, "mc_startMomentum"
		, "mc_truthlength"
		, "track_mcPDG"
		, "n_showers"
		, "sh_start_X"
		, "sh_start_Y"
		, "sh_start_Z"
		, "sh_direction_X"
		, "sh_direction_Y"
		, "sh_direction_Z"
		, "sh_length"
		});
    AllowBranches(allowed);

    // next(mEntry);
}

void EventDisplay::PrintInfo()
{
    mTree->GetEntry(mEntry);
    cout<<"Entry "<<mEntry
	<<", Subrun No: "<<mAnaSel->subRunNo
	<<", Event No: "<<mAnaSel->eventNo<<endl;
    cout<<"Kaon track: "<<endl
	<<"    KinE: "<<(mAnaSel->mc_startMomentum[1][3]-0.4937)*1e3<<" MeV"<<endl
	<<"    Length: "<<mAnaSel->mc_truthlength[1]<<" cm"<<endl;
    cout<<"N Reco Tracks: "<< mAnaSel->n_reco_tracks<<endl;

    // TPC with hits:

}

int EventDisplay::GetEntry(int entry)
{
    mEntry = entry;
    return mTree->GetEntry(mEntry);
}

int EventDisplay::next(const int entry)
{
    if (entry == -1)
	mEntry++;
    else
	mEntry = entry;

    PrintInfo();
    plot();

    return 0;
}

int EventDisplay::plot()
{
    // plot mc truth tracks
    static TMultiGraph** mg_mc = 0;
    static TMultiGraph** mg_mc_wv = 0;
    // plot calorimetric points by tracks
    static TMultiGraph** mg_rec = 0;
    static TMultiGraph** mg_rec_wv = 0;

    //mEntry = mTree->Get

    if (mg_mc)
	delete[] mg_mc;
    if (mg_mc_wv)
	delete[] mg_mc_wv;
    if (mg_rec)
	delete[] mg_rec;
    if (mg_rec_wv)
	delete[] mg_rec_wv;

    mg_mc = plotMCparts();
    mg_mc_wv = plotMCpartsWV();
    mg_rec = plotTracks();
    mg_rec_wv = plotTracksWV();

    vector<vector<TArrow*> > showers = GetShowers();

    for (int i = 0; i < 3; i++) {
    	mCanvas -> cd(i+2);
	gPad->Clear();
	gPad->Update();

    	mg_mc[i] -> Draw("APL PMC PLC");

    	mg_rec[i] -> Draw("P PMC PLC");

	gPad->Modified();
	gPad->Update();

	// mg_mc[i] -> Draw("APL");
	// mg_rec[i] -> Draw("P");
	mg_mc[i] -> Draw("PL");
	gPad->Update();

	mCanvasWv->cd(i+1);
	gPad->Clear();
	gPad->Update();
	mg_mc_wv[i] -> Draw("A PL PMC PLC");
	mg_rec_wv[i] -> Draw("P PMC PLC");
	mg_mc_wv[i] -> Draw("PL");

	// // test drift-time axis and correct the graphical axis if necessary
	// if ((mg_mc[0]->GetXaxis()->GetXmin() + mg_mc[0]->GetXaxis()->GetXmax()) < 0) {
	//     // most of the plot is on the negative side of APA
	//     reverseYaxix(mg_mc_wv[i]);
	// }
	gPad->Update();
    }

    // for (auto gr: *(mg_mc[0]->GetListOfGraphs()))
    // 	cout<<"==== "<<((TGraph*)gr)->GetLineColor()<<endl;

    // draw legend
    int ngraphs = mg_mc[0]->GetListOfGraphs()->GetSize();
    int ncolumns = (ngraphs+4)/5;
    int nrows = (ngraphs<5)?ngraphs:(ngraphs/ncolumns + ngraphs%ncolumns);
    TLegend* legmc = new TLegend(0.98-ncolumns*0.22,1.-nrows*0.06-0.06,0.98,1., "Truth");
    for ( auto gr: *(mg_mc[0]->GetListOfGraphs()) ) legmc->AddEntry(gr, gr->GetTitle(), "lp");
    legmc->SetNColumns(ncolumns);
    legmc->SetFillStyle(0);
    //legmc->SetMargin(0);

	// cout<<"Plot MC legend: "<<endl
	//     <<"  ngraphs: "<<ngraphs<<endl
	//     <<"  ncolumns: "<<ncolumns<<endl
	//     <<"  nrows: "<<nrows<<endl;

    if ( mg_rec[0]->GetListOfGraphs() )
	ngraphs = mg_rec[0]->GetListOfGraphs()->GetSize();
    else
	ngraphs = 0;
    ncolumns = (ngraphs+4)/5;
    nrows = (ngraphs<5)?ngraphs:(ngraphs/ncolumns + ngraphs%ncolumns);
    TLegend* legreco = new TLegend(0.98-ncolumns*0.22,0.5-nrows*0.06-0.06,0.98,0.5, "Reco");
    if (ngraphs)
	for ( auto gr: *(mg_rec[0]->GetListOfGraphs()) ) legreco->AddEntry(gr, gr->GetTitle(), "lp");
    legreco->SetNColumns(ncolumns);
    legreco->SetFillStyle(0);

    cout<<"Plot Reco legend: "<<endl
        <<"  ngraphs: "<<ngraphs<<endl
        <<"  ncolumns: "<<ncolumns<<endl
        <<"  nrows: "<<nrows<<endl;

    mCanvas->cd(1);
    gPad->Clear();
    legmc->Draw();
    legreco->Draw();


    // correct colors of 3D view
    int igr = 0;
    for (auto primitive: *(mCanvas3d->GetListOfPrimitives()) ) {
	if ( !strstr(primitive->ClassName(), "TGraph2D") ) continue;

	auto gr = (TGraph*)mg_mc[0]->GetListOfGraphs()->At(igr);

	((TGraph2D*)primitive)->SetLineColor(gr->GetLineColor());
	((TGraph2D*)primitive)->SetMarkerColor(gr->GetMarkerColor());

	// cout<<"==== "<<((TGraph2D*)primitive)->ClassName()<<": "
	//     <<((TGraph2D*)primitive)->GetLineColor()<<endl;

	igr++;
    }
    mCanvas3d->Modified(); mCanvas3d->Update();
    // 	cout<<primitive->GetName()<<endl;



    // draw the showers
    int iview = 0;
    for (auto view: showers) {
	iview++;
	mCanvas->cd(iview);
	for (auto shower: view) {
	    shower->Draw();
	}
    }

    if (showers.size() > 0)
	cout<<"Have "<< showers[0].size() << " reconstructed showers."<<endl;
    else
	cout<<"Have 0 reconstructed showers."<<endl;

    // TPC Info
    double xmin = mg_mc[1] -> GetXaxis()->GetXmin();
    double xmax = mg_mc[1] -> GetXaxis()->GetXmax();
    double ymin = mg_mc[2] -> GetYaxis()->GetXmin();
    double ymax = mg_mc[2] -> GetYaxis()->GetXmax();
    double zmin = mg_mc[2] -> GetXaxis()->GetXmin();
    double zmax = mg_mc[2] -> GetXaxis()->GetXmax();

    double mean_x = (xmin + xmax)*0.5;
    double mean_y = (ymin + ymax)*0.5;
    double mean_z = (zmin + zmax)*0.5;

    static const double tpc_width = 232.4;
    static const double tpc_offset = -0.9;

    int tpc = (mean_z - tpc_offset)/tpc_width;
    tpc *= 4;
    tpc += (mean_x > 0);
    tpc += (mean_y > 0) * 2;

    cout<<"Event centered in TPC "<<tpc<<endl;

    mCanvas->Update();
    return 0;
}

void EventDisplay::AllowBranches(const std::vector<string>& br_list)
{
    for (auto it = br_list.begin(); it != br_list.end(); it++) {
	mTree -> SetBranchStatus(it->c_str(),1);
    }
    return;
}

TMultiGraph** EventDisplay::plotMCparts()
{
    static TDatabasePDG* pdgdb = TDatabasePDG::Instance();

    int nMCparts = mAnaSel->mc_npart;

    TMultiGraph **mg = new TMultiGraph*[3];
    vector<TGraph2D*> grv;

    const char* titles[3] = {";X [cm];Y [cm]", ";X [cm];Z [cm]", ";Z [cm];Y [cm]"};


    // // find kaon track
    // int k_id = -1;
    // int k_idx = -1;
    // for (int itrk = 0; itrk < nMCparts; itrk++) {
    // 	if ( mAnaSel -> mc_pdg[itrk] == 321 ) {
    // 	    k_idx = itrk;
    // 	    k_id = mAnaSel -> mc_id[itrk];
    // 	}
    // }

    // multigraphs for 2D projections and wire-time views
    for (int i = 0; i < 3; i++) {
	mg[i] = new TMultiGraph(Form("mgMCparts%d",i), Form("Truth Particles%s", titles[i]));
    }

    cout<<"Going through "<<nMCparts<<" MC particles"<<endl;

    int icolor = 2;
    for (int itrk = 0; itrk < nMCparts; itrk++) {
	// don't show neutrinos, neutrons, heavy ions, gammas
	if (abs(mAnaSel->mc_pdg[itrk]) == 12 ||
	    abs(mAnaSel->mc_pdg[itrk]) == 14 ||
	    abs(mAnaSel->mc_pdg[itrk]) == 16 ||
	    mAnaSel->mc_pdg[itrk] == 2112 ||
	    mAnaSel->mc_pdg[itrk] > 100000 ||
	    mAnaSel->mc_pdg[itrk] == 22 ) continue;
	auto gr2d = new TGraph2D();
	grv.push_back(gr2d);

	vector<double> x[3], y[3]; // projections xy, xz, zy
	int npoints = mAnaSel->mc_npoints[itrk];
	for (int ipt = 0; ipt < npoints; ipt++) {
	    double* coord = (double*)mAnaSel->mc_xyz[itrk][ipt];
	    x[0].push_back(coord[0]);
	    y[0].push_back(coord[1]);
	    x[1].push_back(coord[0]);
	    y[1].push_back(coord[2]);
	    x[2].push_back(coord[2]);
	    y[2].push_back(coord[1]);
	    // do 3d pragphs; put z, x, y as graph's x, y, z, respectively
	    gr2d->SetPoint(ipt, coord[2], coord[0], coord[1]);
	    gr2d->SetMarkerColor(icolor);
	    gr2d->SetLineColor(icolor);
	}

	// create and plot the 3 graphs
	int pdg = mAnaSel->mc_pdg[itrk];
	TString name = pdgdb->GetParticle(pdg)->GetName();
	if (name.Contains("gamma"))
	    name.ReplaceAll("gamma", "#gamma");
	for (int i = 0; i < 3; i++) {
	    auto gr = new TGraph(x[i].size(), x[i].data(), y[i].data());
	    gr->SetName(Form("grMC%d_p%d", itrk, i));
	    gr->SetTitle(name);
	    //gr->SetTitle(Form("%d", pdg));
	    gr->SetLineWidth(3);
	    gr->SetMarkerSize(1.);
	    gr->SetMarkerColor(icolor);
	    gr->SetMarkerStyle(kFullDiamond);
	    gr->SetLineColor(icolor);
	    mg[i]->Add(gr);
	}
	icolor++;
    }

    cout<<"Found and to be plotted "<<grv.size()<<" particles"<<endl;

    plot2Dgraphs(grv);

    return mg;
}

TMultiGraph** EventDisplay::plotMCpartsWV()
{
    static TDatabasePDG* pdgdb = TDatabasePDG::Instance();

    int nMCparts = mAnaSel->mc_npart;

    TMultiGraph **mg = new TMultiGraph*[3];

    const char* titles[3] = {"Collection Plane;Z direction [cm];X [cm]",
			     "U Plane;U direction [cm];X [cm]",
			     "V Plane;V direction [cm];X [cm]"};

    // multigraphs for 2D projections and wire-time views
    for (int i = 0; i < 3; i++) {
	mg[i] = new TMultiGraph(Form("mgMCparts%d",i), Form("%s", titles[i]));
    }

    cout<<"Going through "<<nMCparts<<" MC particles"<<endl;

    int icolor = 2;
    for (int itrk = 0; itrk < nMCparts; itrk++) {
	// don't show neutrinos, neutrons, heavy ions, gammas
	if (abs(mAnaSel->mc_pdg[itrk]) == 12 ||
	    abs(mAnaSel->mc_pdg[itrk]) == 14 ||
	    abs(mAnaSel->mc_pdg[itrk]) == 16 ||
	    mAnaSel->mc_pdg[itrk] == 2112 ||
	    mAnaSel->mc_pdg[itrk] > 100000 ||
	    mAnaSel->mc_pdg[itrk] == 22 ) continue;

	vector<double> x[3], y[3]; // projections xy, xz, zy
	int npoints = mAnaSel->mc_npoints[itrk];
	for (int ipt = 0; ipt < npoints; ipt++) {
	    double* coord = (double*)mAnaSel->mc_xyz[itrk][ipt];
	    TVector3 space_point(coord);
	    double u = space_point.Dot(u_dir);
	    double v = space_point.Dot(v_dir);
	    double z = coord[2];
	    double drift = coord[0];
	    // if (coord[0] > 0) {
	    // 	double tmp = u;
	    // 	u = v;
	    // 	v = tmp;
	    // } else
	    // 	drift = -drift;

	    x[0].push_back(coord[2]);
	    y[0].push_back(drift);
	    x[1].push_back(u);
	    y[1].push_back(drift);
	    x[2].push_back(v);
	    y[2].push_back(drift);
	}

	// create and plot the 3 graphs
	int pdg = mAnaSel->mc_pdg[itrk];
	TString name = pdgdb->GetParticle(pdg)->GetName();
	if (name.Contains("gamma"))
	    name.ReplaceAll("gamma", "#gamma");
	for (int i = 0; i < 3; i++) {
	    auto gr = new TGraph(x[i].size(), x[i].data(), y[i].data());
	    gr->SetName(Form("grMCwv%d_p%d", itrk, i));
	    gr->SetTitle(name);
	    //gr->SetTitle(Form("%d", pdg));
	    gr->SetLineWidth(3);
	    gr->SetMarkerSize(1.);
	    gr->SetMarkerColor(icolor);
	    gr->SetMarkerStyle(kFullDiamond);
	    gr->SetLineColor(icolor);
	    mg[i]->Add(gr);
	}
	icolor++;
    }

    return mg;
}

TH1* EventDisplay::plot2Dgraphs(vector<TGraph2D*>& grv)
{
    mCanvas3d -> cd();

    vector<double> xmins, ymins, zmins, xmaxs, ymaxs, zmaxs;
    int size = grv.size();
    if (size == 0) return 0;
    for(int i = 0; i < size; i++) {
	xmins.push_back(grv[i]->GetXmin());
	ymins.push_back(grv[i]->GetYmin());
	zmins.push_back(grv[i]->GetZmin());
	xmaxs.push_back(grv[i]->GetXmax());
	ymaxs.push_back(grv[i]->GetYmax());
	zmaxs.push_back(grv[i]->GetZmax());

    }

    double xmin = TMath::MinElement(xmins.size(), xmins.data());
    double xmax = TMath::MaxElement(xmaxs.size(), xmaxs.data());
    double ymin = TMath::MinElement(ymins.size(), ymins.data());
    double ymax = TMath::MaxElement(ymaxs.size(), ymaxs.data());
    double zmin = TMath::MinElement(zmins.size(), zmins.data());
    double zmax = TMath::MaxElement(zmaxs.size(), zmaxs.data());

    // cout<<xmin<<", "<<xmax<<endl
    // 	<<ymin<<", "<<ymax<<endl
    // 	<<zmin<<", "<<zmax<<endl;

    TH1* hold = (TH1*)gDirectory->Get("hgr2d");
    if (hold) delete hold;

    auto hnew = new TH3F("hgr2d","", 1, xmin, xmax, 1, ymin, ymax, 1, zmin, zmax);
    hnew->SetStats(0);
    hnew->Draw();
    for(int i = 0; i < size; i++) {
	    grv[i]->Draw("p line same pmc plc");
    }

    hnew->SetTitle(";z [cm];x [cm]; y [cm]");
    hnew->GetXaxis()->SetTitleOffset(1.2);
    hnew->GetYaxis()->SetTitleOffset(1.2);
    hnew->GetZaxis()->SetTitleOffset(1.2);

    //cout<<grv[0]->GetHistogram()->GetName()<<endl;

    //gPad->BuildLegend();

    return hnew;
}

TMultiGraph** EventDisplay::plotTracks() {
    static TDatabasePDG* pdgdb = TDatabasePDG::Instance();

    int nTracks  = mAnaSel->n_reco_tracks;

    TMultiGraph **mg = new TMultiGraph*[3];

    const char* titles[3] = {" - Front;X [cm];Y [cm]", "- Top;X [cm];Z [cm]", " - Side;Z [cm];Y [cm]"};

    for (int i = 0; i < 3; i++) {
	mg[i] = new TMultiGraph(Form("mgTracks%d",i), Form("Reco Tracks%s", titles[i]));
    }

    for (int itrk = 0; itrk < nTracks; itrk++) {
	vector<double> x[3], y[3]; // projections xy, xz, zy
	// for (int ipln = 0; ipln < 3; ipln++) {
	double* reco_pts_xyz = (double*)mAnaSel->track_point_xyz[itrk];
	int npoints = mAnaSel->n_track_points[itrk];
	for (int ipts = 0; ipts < npoints; ipts++) {
	    x[0].push_back(reco_pts_xyz[ipts*3]);
	    y[0].push_back(reco_pts_xyz[ipts*3 + 1]);
	    x[1].push_back(reco_pts_xyz[ipts*3]);
	    y[1].push_back(reco_pts_xyz[ipts*3 + 2]);
	    x[2].push_back(reco_pts_xyz[ipts*3 + 2]);
	    y[2].push_back(reco_pts_xyz[ipts*3 + 1]);
	}


	// create and plot the 3 graphs
	static const int viewX[3] = {0, 0, 2};
	static const int viewY[3] = {1, 2, 1};
	double* vertex = mAnaSel->track_vtx[itrk];
	double* end = mAnaSel->track_end[itrk];
	int pdg = mAnaSel->track_mcPDG[itrk];
	TString name = pdgdb->GetParticle(pdg)->GetName();
	// cout<<"Track "<<itrk<<" vertex: "
	//     <<vertex[0]<<", "<<vertex[1]<<", "<<vertex[2]<<", "
	//     <<"end: "
	//     <<end[0]<<", "<<end[1]<<", "<<end[2]<<endl;
	if (name.Contains("gamma"))
	    name.ReplaceAll("gamma", "#gamma");
	for (int i = 0; i < 3; i++) {
	    auto gr = new TGraph(x[i].size(), x[i].data(), y[i].data());
	    gr->SetName(Form("grTrack%d_p%d", itrk, i));
	    gr->SetTitle(name);
	    //gr->SetTitle(Form("%d", pdg));
	    gr->SetMarkerSize(1.);
	    gr->SetMarkerStyle(kOpenSquare);
	    gr->SetMarkerColor(itrk+2);
	    gr->SetLineColor(itrk+2);

	    // add tracks direction
	    auto arrow = new TArrow(vertex[viewX[i]], vertex[viewY[i]],
				    end[viewX[i]], end[viewY[i]], 0.015, "|-|>");
	    arrow->SetNDC(kFALSE);
	    gr->GetListOfFunctions()->Add(arrow);

	    mg[i]->Add(gr, "p");
	}
    }

    return mg;
}

TMultiGraph** EventDisplay::plotTracksWV() {
    static TDatabasePDG* pdgdb = TDatabasePDG::Instance();

    int nTracks  = mAnaSel->n_reco_tracks;

    TMultiGraph **mg = new TMultiGraph*[3];

    //const char* titles[3] = {" - Front;X [cm];Y [cm]", "- Top;X [cm];Z [cm]", " - Side;Z [cm];Y [cm]"};
    const char* titles[3] = {";Z direction [cm];X [cm]", ";U direction [cm];X [cm]", ";V direction [cm];X [cm]"};

    for (int i = 0; i < 3; i++) {
	mg[i] = new TMultiGraph(Form("mgTracksWv%d",i), Form("Reco Tracks%s", titles[i]));
    }

    for (int itrk = 0; itrk < nTracks; itrk++) {
	vector<double> x[3], y[3]; // projections xy, xz, zy
	// for (int ipln = 0; ipln < 3; ipln++) {
	double* reco_pts_xyz = (double*)mAnaSel->track_point_xyz[itrk];
	int npoints = mAnaSel->n_track_points[itrk];
	for (int ipts = 0; ipts < npoints; ipts++) {
	    double* coord = reco_pts_xyz + ipts*3;
	    TVector3 space_point(coord);
	    double u = space_point.Dot(u_dir);
	    double v = space_point.Dot(v_dir);
	    double z = coord[2];
	    double drift = coord[0];
	    // if (coord[0] > 0) {
	    // 	double tmp = u;
	    // 	u = v;
	    // 	v = tmp;
	    // } else
	    // 	drift = -drift;

	    x[0].push_back(coord[2]);
	    y[0].push_back(drift);
	    x[1].push_back(u);
	    y[1].push_back(drift);
	    x[2].push_back(v);
	    y[2].push_back(drift);
	}


	// create and plot the 3 graphs
	static const int viewX[3] = {0, 0, 2};
	static const int viewY[3] = {1, 2, 1};
	double* vertex = mAnaSel->track_vtx[itrk];
	double* end = mAnaSel->track_end[itrk];
	int pdg = mAnaSel->track_mcPDG[itrk];
	TString name = pdgdb->GetParticle(pdg)->GetName();
	// cout<<"Track "<<itrk<<" vertex: "
	//     <<vertex[0]<<", "<<vertex[1]<<", "<<vertex[2]<<", "
	//     <<"end: "
	//     <<end[0]<<", "<<end[1]<<", "<<end[2]<<endl;
	if (name.Contains("gamma"))
	    name.ReplaceAll("gamma", "#gamma");
	for (int i = 0; i < 3; i++) {
	    auto gr = new TGraph(x[i].size(), x[i].data(), y[i].data());
	    gr->SetName(Form("grTrack%d_p%d", itrk, i));
	    gr->SetTitle(name);
	    //gr->SetTitle(Form("%d", pdg));
	    gr->SetMarkerSize(1.);
	    gr->SetMarkerStyle(kOpenSquare);
	    gr->SetMarkerColor(itrk+2);
	    gr->SetLineColor(itrk+2);

	    // // add tracks direction
	    // auto arrow = new TArrow(vertex[viewX[i]], vertex[viewY[i]],
	    // 			    end[viewX[i]], end[viewY[i]], 0.015, "|-|>");
	    // arrow->SetNDC(kFALSE);
	    // gr->GetListOfFunctions()->Add(arrow);

	    mg[i]->Add(gr, "p");
	}
    }

    return mg;
}

std::vector<std::vector<TArrow*> > EventDisplay::GetShowers()
{
    std::vector<std::vector<TArrow*> > showers;
    showers.resize(3);

    int size = mAnaSel->n_showers;
    if (size >  MAX_SHOWERS)
	size = MAX_SHOWERS;

    //cout<<"Will go through "<<size<<" showers"<<endl;

    for (int ishwr = 0; ishwr < size ; ishwr++) {
	TVector3 shower(mAnaSel->sh_direction_X[ishwr],
			mAnaSel->sh_direction_Y[ishwr],
			mAnaSel->sh_direction_Z[ishwr]);
	shower.SetMag(mAnaSel->sh_length[ishwr]);
	TVector3 start(mAnaSel->sh_start_X[ishwr],
		       mAnaSel->sh_start_Y[ishwr],
		       mAnaSel->sh_start_Z[ishwr]);
	TVector3 end(shower+start);

	// create 3 views, xy, xz, zy
	double startX[3] = {start.X(), start.X(), start.Z()};
	double startY[3] = {start.Y(), start.Z(), start.Y()};
	double endX[3] = {end.X(), end.X(), end.Z()};
	double endY[3] = {end.Y(), end.Z(), end.Y()};

	//cout<<"Will fill views for shower "<<ishwr<<endl;

	for (int iview = 0; iview < 3; iview++) {
	    auto arrow = new TArrow(startX[iview], startY[iview],
				    endX[iview], endY[iview], 0.015, "|-|>");
	    arrow->SetLineWidth(2);
	    showers[iview].push_back(arrow);
	}
    }

    return showers;
}

void EventDisplay::reverseYaxix(TMultiGraph* mg)
{
    TAxis* ax_orig = mg->GetYaxis();
    ax_orig->SetLabelOffset(999);
    ax_orig->SetTickLength(0);

    gPad->Update();
    auto ax_new = new TGaxis(gPad->GetUxmin(),
			     gPad->GetUymax(),
			     gPad->GetUxmin()-0.001,
			     gPad->GetUymin(),
			     -ax_orig->GetXmax(),
			     -ax_orig->GetXmin(),
			     510,"+");
    ax_new->SetLabelOffset(-0.03);
    ax_new->Draw();
}


void NDKAna::Loop()
{}

#endif
