{
  TCanvas* canvas = new TCanvas("canvas");

  TFile* infile1 = new TFile("/dune/app/users/tstokes/v09_22_02/srcs/protoduneana/protoduneana/PDK/ndk_test.root");
  infile1->cd("ndkTest");
  infile1->ls(); // prints out what histogram names are in there

  gStyle->SetOptStat(0);
  //TH2F* lengthDiff_HIST = (TH2F*) infile1->Get("lsuanalyzer/lengthDiff_GRAPH");
  TH1F* mc_PDG = (TH1F*) infile1->("ndkTest/mc_PDG");
  
  lengthDiff_HIST->RebinY(10);
  lengthDiff_HIST->Draw("colorz");

  canvas->SaveAs("lengthDiff_HIST.png");
}

//root -b -q -x plotter.C
//This is just an example of how to work the plotter.C macro.  Change to fit your needs
