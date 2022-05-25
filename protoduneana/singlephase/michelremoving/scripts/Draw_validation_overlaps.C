void Draw_validation_overlaps(int plane){
  
  TString histname = "dedx_X";
  //TString histname = "dqdx_X_correction";


  TString plane_str = Form("%d", plane);
  cout << "Plane " << plane_str << endl;
  map<TString, TH1F*> maphist;
  map<TString, TFile*> mapfile;

  for(int i = 0; i < 11; i++){
    double ref_dedx = 1.5 + 0.1 * i;
    TString ref_dedx_str = Form("%g", ref_dedx);
    if(i == 5) ref_dedx_str = ref_dedx_str + ".0";
    cout << ref_dedx_str << endl;
    
    mapfile[ref_dedx_str] = new TFile("./outputs/validation_mich2_run_valid_" + ref_dedx_str + ".root");
    //mapfile[ref_dedx_str] = new TFile("./outputs/Xcalo_mich2_r5387_" + ref_dedx_str + ".root");
    maphist[ref_dedx_str] = (TH1F*) gDirectory -> Get(histname + "_hist_" + plane_str);
  }

  /*
  mapfile["BOX"]  = new TFile("./outputs/Xcalo_mich2_r5387.root");
  maphist["BOX"] = (TH1F*) gDirectory -> Get(histname + "_hist_" + plane_str);
  */

  TCanvas *c = new TCanvas("", "", 800, 600);
  TH1F* h_template = new TH1F("", "", 144, -360, 360);
  if (histname.Contains("dedx")) h_template -> GetYaxis() -> SetRangeUser(1.86, 1.94);
  if (histname.Contains("dqdx")) h_template -> GetYaxis() -> SetRangeUser(56500, 60000);
  if (histname.Contains("correction")) h_template -> GetYaxis() -> SetRangeUser(0.9, 1.1);
  h_template -> SetLineColor(kBlack);
  h_template -> SetTitle("Plane " + plane_str);
  h_template -> SetLineWidth(2);
  h_template -> SetStats(0);
  if (histname.Contains("dedx"))h_template -> GetYaxis() -> SetTitle("dE/dx (MeV/cm)");
  if (histname.Contains("dqdx"))h_template -> GetYaxis() -> SetTitle("dQ/dx (e/cm)");
  h_template -> GetXaxis() -> SetTitle("X coordinate (cm)");
  h_template -> Draw();

  TLegend *l = new TLegend(0.40, 0.15, 0.60, 0.40);
  l-> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(1001);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);

  Int_t colour_array[] = {600, 591, 616, 611, 632, 625, 400, 394, 416, 418, 432, 426};

  for(int i = 0; i < 11; i++){
    //if(i == 5) continue;
    double ref_dedx = 1.5 + 0.1 * i;
    TString ref_dedx_str = Form("%g", ref_dedx);
    if(i == 5) ref_dedx_str = ref_dedx_str + ".0";
    maphist[ref_dedx_str] -> SetLineColor(colour_array[i]);
    maphist[ref_dedx_str] -> SetLineWidth(2);
    maphist[ref_dedx_str] -> Draw("same");
    l -> AddEntry(maphist[ref_dedx_str], "Ref. dE/dx = " + ref_dedx_str + " MeV/cm", "l");
  }

  /*
  maphist["BOX"] -> SetLineColor(kBlack);
  maphist["BOX"] -> SetLineWidth(2);
  maphist["BOX"] -> Draw("same");
  l -> AddEntry(maphist["BOX"], "Box Ref. dE/dx = 1.9 MeV/cm", "l");
  */

  l -> Draw("same");

  c -> SaveAs("./plots/" + histname + "_Birks_plane" + plane_str + ".pdf");

  
  for(int i = 0; i < 11; i++){
    double ref_dedx = 1.5 + 0.1 * i;
    TString ref_dedx_str = Form("%g", ref_dedx);
    if(i == 5) ref_dedx_str = ref_dedx_str + ".0";
    mapfile[ref_dedx_str] -> Close();
  }
  //mapfile["BOX"] -> Close();
}
