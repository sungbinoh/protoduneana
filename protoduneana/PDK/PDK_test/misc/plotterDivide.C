void {
  TCanvas* canvas = new TCanvas("canvas","canvas",900,800);
  //canvas->Divide(2,1);

  TFile* infile1 = new TFile("/dune/app/users/hmeyer5/v07_13_00/srcs/dunetpc/dune/LSU/LSU_hist.root");
  infile1->cd("lsuanalyzer");
  //infile1->ls(); // prints out what histogram names are in there
  
  TH2F* hist = (TH2F*) infile1->Get("lsuanalyzer/recoPolyFractionalDiffvsTrueMomentum_HIST");
  //TH2F* thetaYZprimevsPolygonalSegmentNumber_HIST = (TH2F*) infile1->Get("lsuanalyzer/truePolyThetaYZprimevsSegmentNumber_HIST");
  //canvas->cd(1);
  TProfile* profile1 = hist->ProfileX("Theta XZ prime vs Polygonal Segment Number Profile",1,-1,"");
  //hist->Draw("colz");
  profile1->Rebin(10);
  profile1->SetLineColor(kRed);
  profile1->SetLineWidth(1);
  profile1->GetXaxis()->SetTitle("True Momentum (GeV)");
  profile1->GetYaxis()->SetTitle("Fractional Bias");
  profile1->SetTitle("Reco Poly Fractional Inverse Momentum Bias vs. True Momentum");
  //profile1->Draw("same");

  TH1F* he = new TH1F("recoLinearFractionalDiffvsTrueMomentum_HIST","Reco Linear Inverse Fractional Resolution;True Momentum;Inverse Fractional Difference",hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  Int_t n = hist->GetNcells();
  for(Int_t i = 0; i < n; i++)
    {
      //he->SetBinContent(i*10,(profile1->GetBinError(i))/sqrt(n));
      he->SetBinContent(i*10,(profile1->GetBinError(i)));
      he->SetBinError(i*10,100/*(profile1->GetBinError(i))/sqrt(2*n)*/);
    }
  canvas->Divide(1,2);
  canvas->cd(2);
  he->SetTitle("Reco Poly Fractional Inverse Momentum Resolution");
  he->GetXaxis()->SetTitle("True Momentum (GeV)");
  he->GetYaxis()->SetTitle("Fractional Resolution");
  he->Draw("HIST *P");
  canvas->cd(1);
  profile1->Draw();

  /*canvas->cd(2);
  TProfile* profile2 = thetaYZprimevsPolygonalSegmentNumber_HIST->ProfileX("Theta YZ prime vs Polygonal Segment Number Profile",1,-1,"");
  profile2->Draw();
  
  /*TH2F* polygonalMCSvsTrueMomentum_HIST = (TH2F*) infile1->Get("lsuanalyzer/linearMCSvsTrueMomentum_HIST");
    polygonalMCSvsTrueMomentum_HIST->Fit("pol1");*/
  
  //canvas->SaveAs("lengthDiff_HIST.png");
}

//root -q -x plotter.C
//This is just an example of how to work the plotter.C macro.  Change to fit your needs
