void Draw_validation(){

  //TFile *f0 = new TFile("validation_mich2_run5387.root");
  TFile *f0 = new TFile("./outputs/validation_mich2_run_valid_1.9.root");

  TH1F *dedx_x_0 = (TH1F*) gDirectory -> Get("dedx_X_hist_0");
  TH1F *dedx_x_1 = (TH1F*) gDirectory -> Get("dedx_X_hist_1");
  TH1F *dedx_x_2 = (TH1F*) gDirectory -> Get("dedx_X_hist_2");

  TCanvas *c = new TCanvas("", "", 800, 600);
  dedx_x_0 -> GetYaxis() -> SetRangeUser(1.88, 1.95);
  dedx_x_0 -> SetLineColor(kBlack);
  dedx_x_0 -> SetTitle("");
  dedx_x_0 -> SetLineWidth(2);
  dedx_x_0 -> SetStats(0);
  dedx_x_0 -> Draw();

  dedx_x_1 -> SetLineColor(kRed);
  dedx_x_1 -> SetLineWidth(2);
  dedx_x_1 -> Draw("same");

  dedx_x_2 -> SetLineColor(kGreen);
  dedx_x_2 -> SetLineWidth(2);
  dedx_x_2 -> Draw("same");

  TLegend *l = new TLegend(0.12, 0.75, 0.25, 0.88);
  l-> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(1001);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);

  l -> AddEntry(dedx_x_0, "Plane 0", "l");
  l -> AddEntry(dedx_x_1, "Plane 1", "l");
  l -> AddEntry(dedx_x_2, "Plane 2", "l");

  l -> Draw("same");

  c -> SaveAs("./plots/dedx_X_Birks_1.9.pdf");


}
