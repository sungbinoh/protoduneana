#include "canvas_margin.h"

void Draw_Birks_parameters(){
  setTDRStyle();

  double A[3] = {0.803304, 0.793268, 0.823802};
  double A_err[3] = {0.050026, 0.051386, 0.055829};
  double kB[3] = {0.050725, 0.044621, 0.064649};
  double kB_err[3] = {0.024850, 0.025833, 0.027421};
  double A_ICARUS = 0.800;
  double A_ICARUS_err = 0.003;
  double kB_ICARUS = 0.0486;
  double kB_ICARUS_err = 0.0006;

  double y[3] = {1., 2., 3.};
  double y_err[3] = {0., 0., 0.};

  TCanvas *c = new TCanvas("", "", 800, 600);
  TPad *pad1 = new TPad("", "", 0.1, 0, 0.55, 0.95);
  //pad1 -> SetTopMargin( 0.07 );
  pad1 -> SetLeftMargin( 0.15 );
  pad1 -> SetRightMargin( 0.03 );
  pad1 -> Draw();
  pad1 -> cd();
  
  TH1F * pad1_template = new TH1F("", "", 1, 0.7, 0.9);
  TString label[3] = {"Plane0", "Plane1", "Plane2"};
  pad1_template -> GetYaxis() -> SetRangeUser(0., 4.);
  pad1_template -> GetYaxis() -> SetLabelSize(0);
  pad1_template -> GetXaxis() -> SetTitle("A");
  pad1_template -> GetXaxis() -> CenterTitle(true);
  pad1_template -> GetXaxis() -> SetNdivisions(2110);
  pad1_template -> GetXaxis() -> SetLabelSize(0.04);
  pad1_template -> SetStats(0);
  pad1_template -> Draw();

  TBox *pad1_box = new TBox(A_ICARUS - A_ICARUS_err, 0., A_ICARUS + A_ICARUS_err, 4.);
  pad1_box -> SetFillColorAlpha(0,0);
  pad1_box -> Draw("same");

  TLine *pad1_line = new TLine(A_ICARUS, 0., A_ICARUS, 4.);
  pad1_line -> SetLineColor(kRed);
  pad1_line -> Draw("same");

  TGraphErrors *A_gr = new TGraphErrors(3, A, y, A_err, y_err);
  A_gr -> Draw("epsame");

  TLegend *legend1 = new TLegend(0.72, 0.85, 0.97, 0.95);
  legend1 -> AddEntry(pad1_line, "ICARUS", "l");
  legend1 -> SetFillStyle(1001);
  //legend1 -> SetLineColor(kBlack);
  //legend1 -> SetBorderSize(1);
  legend1 -> Draw("same");

  gPad->RedrawAxis();

  c -> cd();

  TPad *pad2 = new TPad("", "", 0.55, 0, 1., 0.95);
  //pad2 -> SetTopMargin( 0.07 );
  pad2 -> SetLeftMargin( 0.03 );
  pad2 -> SetRightMargin( 0.15 );
  pad2 -> Draw();
  pad2 -> cd();

  TH1F * pad2_template = new TH1F("", "", 1, 0.00, 0.10);
  pad2_template -> GetYaxis() -> SetRangeUser(0., 4.);
  pad2_template-> GetYaxis() -> SetLabelSize(0);
  pad2_template -> GetXaxis() -> SetTitle("kB");
  pad2_template -> GetXaxis() -> CenterTitle(true);
  pad2_template -> GetXaxis() -> SetNdivisions(1110);
  pad2_template -> GetXaxis() -> SetLabelSize(0.04);
  pad2_template -> SetStats(0);
  pad2_template ->Draw();

  TBox *pad2_box = new TBox(kB_ICARUS - kB_ICARUS_err, 0., kB_ICARUS + kB_ICARUS_err, 4.);
  pad2_box -> SetFillColorAlpha(0,0);
  pad2_box -> Draw("same");

  TLine *pad2_line = new TLine(kB_ICARUS, 0., kB_ICARUS, 4.);
  pad2_line -> SetLineColor(kRed);
  pad2_line -> Draw("same");

  TGraphErrors *kB_gr = new TGraphErrors(3, kB, y, kB_err, y_err);
  kB_gr -> Draw("epsame");

  TLegend *legend2 = new TLegend(0.60, 0.85, 0.85, 0.95);
  legend2 -> AddEntry(pad2_line, "ICARUS", "l");
  legend2 -> SetFillStyle(1001);
  legend2 -> Draw("same");

  gPad->RedrawAxis();

  c -> cd();
  
  TLatex latex_plane0, latex_plane1, latex_plane2;
  double text_size = 0.035;
  latex_plane0.SetNDC();
  latex_plane1.SetNDC();
  latex_plane2.SetNDC();
  latex_plane0.SetTextSize(text_size);
  latex_plane1.SetTextSize(text_size);
  latex_plane2.SetTextSize(text_size);
  latex_plane0.DrawLatex(0.05, 0.30, "Plane0");
  latex_plane0.DrawLatex(0.05, 0.50, "Plane1");
  latex_plane0.DrawLatex(0.05, 0.70, "Plane2");

  TLatex latex_data, latex_particle;
  latex_data.SetNDC();
  latex_particle.SetNDC();
  latex_data.SetTextSize(0.04);
  latex_particle.SetTextSize(0.035);
  latex_data.DrawLatex(0.17, 0.92, "ProtoDUNE-SP #font[42]{#it{Run 5387}}");
  latex_particle.DrawLatex(0.57, 0.92, "#font[42]{Stopping} #mu #font[42]{(#it{KE}} #font[42]{= [50, 450] #it{MeV})}");
  

  c -> SaveAs("./plots/Recom_comparison_Birks.pdf");

}
