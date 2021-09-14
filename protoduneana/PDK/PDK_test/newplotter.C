{
  TCanvas* canvas = new TCanvas("canvas");
  TFile* infile1 = new TFile("/dune/app/users/tstokes/v09_22_02/srcs/protoduneana/protoduneana/PDK/ndk_test.root");
  infile1->cd("NDKAna");
  infile1->ls(); // prints out what histogram names are in there

  gStyle->SetOptStat(0);
  //True Histograms
  TH1F* MC_length         = (TH1F*) infile1->Get("NDKAna/h_MC_length");
  TH1F* MC_length_Kaon   = (TH1F*) infile1->Get("NDKAna/h_MC_length_Kaons");
  TH1F* MC_length_Muon   = (TH1F*) infile1->Get("NDKAna/h_MC_length_Muons");
  TH1F* MC_length_Pion   = (TH1F*) infile1->Get("NDKAna/h_MC_length_Pions");
  //Reco Histgrams
  TH1F* track_length_Muon = (TH1F*) infile1->Get("NDKAna/h_track_length_Muons");
  TH1F* track_length_Kaon = (TH1F*) infile1->Get("NDKAna/h_track_length_Kaons");
  TH1F* track_length_Pion = (TH1F*) infile1->Get("NDKAna/h_track_length_Pions");  
  TH1F* track_length      = (TH1F*) infile1->Get("NDKAna/h_track_length");
 

  track_length_Muon->SetLineColor(kBlue);
  track_length_Muon->Draw();
  //MC_length_Kaon->SetLineColor(kRed);
  //MC_length_Kaon->Draw("same");
  track_length_Muon->Add(MC_length_Pion);
  track_length_Muon->Add(MC_length_Kaon);
  track_length->SetLineColor(kBlack);
  track_length->Draw("same");
  //MC_length_Pion->SetLineColor(kBlack);
  //MC_length_Pion->Draw("same");
  track_length_Kaon->SetTitle("Reco KMP vs Total Reco Tracks");
  track_length_Muon->GetXaxis()->SetTitle("Track Length (cm)");
  track_length_Muon->GetYaxis()->SetTitle("Count");
  track_length_Muon->GetYaxis()->SetRange(0, 55);
  auto legend = new TLegend(0.9, 0.7, 0.7, 0.9);
  legend->SetHeader("Legend", "C");
  legend->AddEntry(track_length_Muon, "KMP", "l");
  legend->AddEntry(track_length, "All Reco", "l");
  //legend->AddEntry(MC_length_Pion, "Pion", "l");  
  legend->Draw();
  canvas->SaveAs("Reco_KMP_vs_All.png");

/*
  track_length_Muon->SetLineColor(kRed);
  track_length_Muon->Add(track_length_Kaon);
  track_length_Muon->Add(track_length_Pion);
  track_length_Muon->Draw();
  track_length->SetLineColor(kBlue);
  track_length->Draw("same");

  track_length_Muon->SetTitle("KMP vs Total Reco Track Lengths");
  track_length_Muon->GetXaxis()->SetTitle("Track Length (cm)");
  track_length_Muon->GetYaxis()->SetTitle("Count");
 
  auto legend = new TLegend(0.9, 0.7, 0.7, 0.9);
  legend->SetHeader("Legend", "C");
  legend->AddEntry(track_length_Muon, "KMP", "l");
  legend->AddEntry(track_length, "Total Reco", "l");
  legend->Draw();
  canvas->SaveAs("KMP_TotalReco.png");
*/
}

//root -b -q -x plotter.C
//This is just an example of how to work the plotter.C macro.  Change to fit your needs
