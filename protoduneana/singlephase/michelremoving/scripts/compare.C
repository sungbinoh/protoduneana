void compare(){

  TFile *f0 = new TFile("Xcalo_mich2_r5387.root");
  TH1F *f0_plane0 = (TH1F*)gDirectory -> Get("dqdx_X_correction_hist_0");
  TH1F *f0_plane1 = (TH1F*)gDirectory -> Get("dqdx_X_correction_hist_1");
  TH1F *f0_plane2 = (TH1F*)gDirectory -> Get("dqdx_X_correction_hist_2");


  TFile *f1 = new TFile("/pnfs/dune/persistent/physicsgroups/protodune/Calibration/Prod4/data/Xcalo_mich2_r5387.root");
  TH1F *f1_plane0 = (TH1F*)gDirectory -> Get("dqdx_X_correction_hist_0");
  TH1F *f1_plane1 = (TH1F*)gDirectory -> Get("dqdx_X_correction_hist_1");
  TH1F *f1_plane2 = (TH1F*)gDirectory -> Get("dqdx_X_correction_hist_2");
  
  TCanvas *c = new TCanvas("", "", 800, 800);
  // == plane 0
  f0_plane0 -> GetYaxis() -> SetRangeUser(0.9, 1.10);
  f0_plane0 -> SetLineColor(kRed);
  f0_plane0 -> SetLineWidth(2);
  f0_plane0 -> SetStats(0);
  f0_plane0 -> Draw();

  f1_plane0 -> SetLineColor(kBlue);
  f1_plane0 -> SetLineWidth(2);
  f1_plane0 -> Draw("same");

  c -> SaveAs("./comparison_0.pdf");

  // == plane 1
  f0_plane1 -> GetYaxis() -> SetRangeUser(0.9, 1.10);
  f0_plane1 -> SetLineColor(kRed);
  f0_plane1 -> SetLineWidth(2);
  f0_plane1 -> SetStats(0);
  f0_plane1 -> Draw();

  f1_plane1 -> SetLineColor(kBlue);
  f1_plane1 -> SetLineWidth(2);
  f1_plane1 -> Draw("same");

  c -> SaveAs("./comparison_1.pdf");

  // == plane 2
  f0_plane2 -> GetYaxis() -> SetRangeUser(0.9, 1.10);
  f0_plane2 -> SetLineColor(kRed);
  f0_plane2 -> SetLineWidth(2);
  f0_plane2 -> SetStats(0);
  f0_plane2 -> Draw();

  f1_plane2 -> SetLineColor(kBlue);
  f1_plane2 -> SetLineWidth(2);
  f1_plane2 -> Draw("same");

  c -> SaveAs("./comparison_2.pdf");

  f0 -> Close();
  f1 -> Close();

}
