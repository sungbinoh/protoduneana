#include "include.h"
#include <iostream>
#include <algorithm>
using namespace std;

void newplotter()
{
  TCanvas* canvas = new TCanvas("canvas");
  TFile* infile1 = new TFile("/dune/app/users/tstokes/PDK_Analysis/srcs/protoduneana/protoduneana/PDK/Output/hA_BR_PDK_Cuts.root");
  infile1->cd("NDKAna");
  infile1->ls(); // prints out what histogram names are in there

  gStyle->SetOptStat(0);
  //True Histograms
  //TH1F* MC_length        =(TH1F*) infile1->Get("NDKAna/h_MC_length");
 // TH1F* MC_length_Kaon   = (TH1F*) infile1->Get("NDKAna/h_MC_length_Kaons");
 // TH1F* MC_length_Muon   = (TH1F*) infile1->Get("NDKAna/h_MC_length_Muons");
  //TH1F* MC_length_Pion   = (TH1F*) infile1->Get("NDKAna/h_MC_length_Pions");
  //TH1F* MC_length_Electron = (TH1F*) infile1->Get("NDKAna/h_MC_length_Electrons");
  //TH1F* MC_length_Proton = (TH1F*) infile1->Get("NDKAna/h_MC_length_Protons");
  //TH1F* MC_length_Neutrino = (TH1F*) infile1->Get("NDKAna/h_MC_length_Neutrinos");

  //Reco Histgrams
  TH1F* track_length_Muon = (TH1F*) infile1->Get("NDKAna/h_track_length_Muons");
  TH1F* track_length_Kaon = (TH1F*) infile1->Get("NDKAna/h_track_length_Kaons");
  //TH1F* track_length_Pion = (TH1F*) infile1->Get("NDKAna/h_track_length_Pions");  
    //TH1F* track_length      = (TH1F*) infile1->Get("NDKAna/h_track_length");
/*
MC_length_Muon->SetLineColor(kBlue);
MC_length_Muon->Add(MC_length_Kaon);
//MC_length_Muon->Add(MC_length_Electron);
//MC_length_Muon->Add(MC_length_Pion);
//MC_length_Muon->Add(MC_length_Proton);
//MC_length_Muon->Add(MC_length_Neutrino);
MC_length_Muon->Draw();
MC_length->SetLineColor(kRed);
MC_length->Draw("same");
MC_length_Muon->SetTitle("Matching PID Tracks to All Tracks");
MC_length_Muon->GetXaxis()->SetTitle("Track Length (cm))");
MC_length_Muon->GetYaxis()->SetTitle("Count");
auto legend = new TLegend(0.9, 0.7, 0.7, 0.9);
legend->SetHeader("Legend", "C");
legend->AddEntry(MC_length_Muon, "PID Tracks", "l");
legend->AddEntry(MC_length, "All True", "l");
//legend->AddEntry(MC_length_Pion, "Pion", "l");
legend->Draw();
canvas->SaveAs("True_PID_vs_All.png");
*/  
 
  track_length_Muon->SetLineColor(kBlue);
  track_length_Muon->Draw();
  track_length_Kaon->SetLineColor(kRed);
  track_length_Kaon->Draw("same");
  //MC_length_Muon->Add(MC_length_Pion);
  //MC_length_Muon->Add(MC_length_Kaon);
  //MC_length_Muon->Draw(); 
  //MC_length->SetLineColor(kBlack);
//  MC_length->Draw("same");
  //MC_length_Pion->SetLineColor(kBlack);
  //MC_length_Pion->Draw("same");
  track_length_Muon->SetTitle("Kaon and Muon Reco Tracks");
  track_length_Muon->GetXaxis()->SetTitle("Track Length (cm)");
  track_length_Muon->GetYaxis()->SetTitle("Count");
  track_length_Muon->GetYaxis()->SetRange(0, 55);
  auto legend = new TLegend(0.9, 0.7, 0.7, 0.9);
  legend->SetHeader("Legend", "C");
  legend->AddEntry(track_length_Muon, "Muon", "l");
  legend->AddEntry(track_length_Kaon, "Kaon", "l");
  //legend->AddEntry(MC_length_Pion, "Pion", "l");  
  legend->Draw();
  canvas->SaveAs("Reco_KaonMuon_2400.png");

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
