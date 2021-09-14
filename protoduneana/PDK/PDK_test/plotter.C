{
  TCanvas* canvas = new TCanvas("canvas");

  TFile* infile1 = new TFile("/dune/app/users/tstokes/v09_22_02/srcs/protoduneana/protoduneana/PDK/Example.root");
  //infile1->cd("");
  infile1->ls(); // prints out what histogram names are in there

 gStyle->SetOptStat(0);
  TH1F* track_length_Muon = (TH1F*) infile1->Get("h_track_length_Muons");
  TH1F* track_length_Kaon = (TH1F*) infile1->Get("h_track_length_Kaons");
  TH1F* track_length_Pion = (TH1F*) infile1->Get("h_track_length_Pions");  
  TH1F* track_length = (TH1F*) infile1->Get("h_track_length");
  
 // TH1F* track_reco_PID = (TH1F*) infile1->Get("h_track_PID_pdg");
 // TH1F* mc_PDG = (TH1F*) infile1->Get("h_mc_pdg");
/*
  mc_PDG->SetLineColor(kRed);
  mc_PDG->Draw();
  track_reco_PID->SetLineColor(kBlue);
  track_reco_PID->Draw("same");
  mc_PDG->SetTitle("True vs Reco PID");
  mc_PDG->GetXaxis()->SetTitle("PID");
  mc_PDG->GetYaxis()->SetTitle("Occurences");
  //canvas->SaveAs("TrueReco_PID.png");
  auto legend = new TLegend(0.9, 0.7, 0.7, 0.9);
  legend->SetHeader("Legend","C");
  legend->AddEntry("mc_PDG", "Truth(red)","l");
  legend->AddEntry("track_reco_PID", "Reco(blue)","l");
  legend->Draw();
  canvas->SaveAs("TrueReco_PID.png");
  */
 

  track_length_Muon->SetLineColor(kRed);
  track_length_Muon->Draw();
  track_length_Muon->Add(track_length_Kaon);  
  track_length_Muon->Add(track_length_Pion); 
  //track_length_Muon->Draw();
  track_length->SetLineColor(kBlue);
  track_length->Draw("same");
//track_lengths->Draw();
//track_length_Kaon->SetLineColor(kGreen);
//track_length_Kaon->Draw("same");
//track_length_Pion->SetLineColor(kBlue);
//track_length_Pion->Draw("same");
track_length_Muon->SetTitle("Reco Track Lengths 2 Tracks");
track_length_Muon->GetXaxis()->SetTitle("Track Length (cm)");
track_length_Muon->GetYaxis()->SetTitle("Occurences");
//  canvas->SaveAs("SHIT.png");

//auto legend = new TLegend(horizontal placement in canvas ,font size? ,width of legend box, 
auto legend = new TLegend(0.9, 0.7, 0.7, 0.9);
legend->SetHeader("Legend","C");
legend->AddEntry(track_length_Muon, "KMP", "l");
legend->AddEntry(track_length, "Total", "l");
//legend->AddEntry(track_length_Pion, "Pions", "l");
legend->Draw();
canvas->SaveAs("SHIT.png");

}

//root -b -q -x plotter.C
//This is just an example of how to work the plotter.C macro.  Change to fit your needs
