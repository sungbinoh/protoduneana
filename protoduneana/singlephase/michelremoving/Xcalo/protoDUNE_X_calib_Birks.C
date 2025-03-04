#define protoDUNE_X_calib_Birks_cxx
#include "protoDUNE_X_calib_Birks.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h> 
#include <string>
#include <TImage.h>
#include <iomanip>
#include <algorithm>
using namespace std;

////defining recombination function
float LAr_density=1.383;
float alp=0.93;
float bet=0.212;
//float dedx=2.08;
float ref_dedx=1.9;
string ref_dedx_str = "1.9";
bool userecom=true;
bool useBirksModel=true;
bool sceon = true;
float recom_factor(float totEf){
  if (!userecom) return 1;
  float xsi=bet*ref_dedx/(LAr_density*totEf);
  float xsi0=bet*ref_dedx/(LAr_density*0.4867);
  float rec0=log(alp+xsi0)/xsi0;
  return (rec0*xsi)/log(alp+xsi);
}
double recom_factor_Birks(float totEF){
  if (!userecom) return 1;
  double Rho = 1.383; // == g/cm^3 (liquid argon density at a pressure 18.0 psia)
  //double A_B = 0.806; //  == ArgoNeuT :  at 0.481kV/cm
  //double k_B = 0.052; // == ArgoNeuT (kV/cm)(g/cm2 )/MeV
  //double A_B = 0.800; // == ICARUS : at 200, 350, 500 V/cm
  double k_B = 0.0486; // == ICARUS (kV/cm)(g/cm2 )/MeV
  
  double dQdx_SCE = 1. / (1 + k_B * ref_dedx / (Rho * totEF));
  double dQdx_nominal = 1. / (1 + k_B * ref_dedx / (Rho * 0.4867));

  //return dQdx_SCE / dQdx_nominal;
  return dQdx_nominal / dQdx_SCE;

}

TFile *ef = TFile::Open("$DUNE_PARDATA_DIR/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
float tot_Ef(float xval,float yval,float zval){
  float E0value=0.4867;
  if (!sceon){
    return E0value;
  }
  else{
    if(xval>=0){
      float ex=E0value+E0value*xpos->GetBinContent(xpos->FindBin(xval,yval,zval));
      float ey=0.0+E0value*ypos->GetBinContent(ypos->FindBin(xval,yval,zval));
      float ez=0.0+E0value*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
      return sqrt(ex*ex+ey*ey+ez*ez);
      // return ex;
    }
    else{
      float ex=E0value+E0value*xneg->GetBinContent(xneg->FindBin(xval,yval,zval));
      float ey=0.0+E0value*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
      float ez=0.0+E0value*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
      return sqrt(ex*ex+ey*ey+ez*ez);
      // return ex;
    }
  }
}

TH3F *zpos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Pos");
TH3F *zneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Neg");
TH3F *zpos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Pos");
TH3F *zneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Neg");

float zoffset(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_fd->GetBinContent(zpos_fd->FindBin(xval,yval,zval));
  }
  else{
  return zneg_fd->GetBinContent(zneg_fd->FindBin(xval,yval,zval));
  }
}
float zoffsetbd(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_bd->GetBinContent(zpos_bd->FindBin(xval,yval,zval));
  }
  else{
  return zneg_bd->GetBinContent(zneg_bd->FindBin(xval,yval,zval));
  }
}




void protoDUNE_X_calib_Birks::Loop(TString mn)
{

  if (fChain == 0) return;

  //int x_bin_size=5;
  //int y_bin_size = 5; // nbiny bins in y direction
  //int z_bin_size = 5; // nbinz bins in z direction
  std::cout<<"efield at the anode neg"<<tot_Ef(-352,300,300)<<std::endl;
  std::cout<<"efield at the anode pos"<<tot_Ef(352,300,300)<<std::endl;

  ///plane_2 details
  TH1F *dqdx_X_hist_2 = new TH1F("dqdx_X_hist_2","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",144,-360,360);
  TH1F *dedx_X_hist_2 = new TH1F("dedx_X_hist_2","plane_2;X Coordinate(cm);dE/dx(MeV/cm)",144,-360,360);
  TH1F *dqdx_X_correction_hist_2 = new TH1F("dqdx_X_correction_hist_2","plane_2;X Coordinate(cm);X Correction factors",144,-360,360);
  TH1F *corrected_dqdx_X_hist_2 = new TH1F("corrected_dqdx_X_hist_2","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",144,-360,360);
  TH1F *dqdx_X_hist_1 = new TH1F("dqdx_X_hist_1","plane_1;X Coordinate(cm);dQ/dx(ADC/cm)",144,-360,360);
  TH1F *dedx_X_hist_1 = new TH1F("dedx_X_hist_1","plane_1;X Coordinate(cm);dE/dx(MeV/cm)",144,-360,360);
  TH1F *dqdx_X_correction_hist_1 = new TH1F("dqdx_X_correction_hist_1","plane_1;X Coordinate(cm);X Correction factors",144,-360,360);
  TH1F *corrected_dqdx_X_hist_1 = new TH1F("corrected_dqdx_X_hist_1","plane_1;X Coordinate(cm);dQ/dx(ADC/cm)",144,-360,360);
  TH1F *dqdx_X_hist_0 = new TH1F("dqdx_X_hist_0","plane_0;X Coordinate(cm);dQ/dx(ADC/cm)",144,-360,360);
  TH1F *dedx_X_hist_0 = new TH1F("dedx_X_hist_0","plane_0;X Coordinate(cm);dE/dx(MeV/cm)",144,-360,360);
  TH1F *dqdx_X_correction_hist_0 = new TH1F("dqdx_X_correction_hist_0","plane_0;X Coordinate(cm);X Correction factors",144,-360,360);
  TH1F *corrected_dqdx_X_hist_0 = new TH1F("corrected_dqdx_X_hist_0","plane_0;X Coordinate(cm);dQ/dx(ADC/cm)",144,-360,360);
  //TH2F *dqdx_vs_X=new TH2F("dqdx_vs_X","dQ/dx vs X for all T0 tagged throughgoing muons;X coordinate(cm);dQ/dx(ADC/cm)",760,-380,380,1000,0,1000);
  //TH1F *max_min=new TH1F("max_min","maximum and minimum values",200,-400,400);
  //TH1F *no_hits=new TH1F("no_hits","no of hits for each bin",200,-400,400);
  TH1F *hdqdx[3];
  TH1F *hdedx[3];
  for (unsigned int i = 0; i<3; ++i){
    hdqdx[i] = new TH1F(Form("hdqdx_%d",i), Form("plane_%d;dQ/dx (ke/cm);Entries",i), 100,0,200);
    hdedx[i] = new TH1F(Form("hdedx_%d",i), Form("plane_%d;dE/dx (MeV/cm);Entries",i), 100,0,10);
  }    

  vector<vector<float>> dqdx_value_2;
  vector<vector<float>> dedx_value_2;
  vector<float> all_dqdx_value_2;
  vector<vector<float>> dqdx_frac_correction_2;
  dqdx_value_2.resize(144);
  dedx_value_2.resize(144);
  dqdx_frac_correction_2.resize(144);
  vector<vector<float>> dqdx_value_1;
  vector<vector<float>> dedx_value_1;
  vector<float> all_dqdx_value_1;
  vector<vector<float>> dqdx_frac_correction_1;
  dqdx_value_1.resize(144);
  dedx_value_1.resize(144);
  dqdx_frac_correction_1.resize(144);
  vector<vector<float>> dqdx_value_0;
  vector<vector<float>> dedx_value_0;
  vector<float> all_dqdx_value_0;
  vector<vector<float>> dqdx_frac_correction_0;
  dqdx_value_0.resize(144);
  dedx_value_0.resize(144);
  dqdx_frac_correction_0.resize(144);

  ////////////////////// Importing Y-Z plane fractional corrections /////////////
  fChain->GetEntry(0);
  if (run>10000) run = 0;
  TFile my_file(Form("YZcalo_mich%s_r%d.root",mn.Data(), run));
  TH2F *YZ_negativeX_hist_2= (TH2F*)my_file.Get("correction_dqdx_ZvsY_negativeX_hist_2");
  TH2F *YZ_positiveX_hist_2= (TH2F*)my_file.Get("correction_dqdx_ZvsY_positiveX_hist_2");
  ////////////////////// Importing Y-Z plane fractional corrections /////////////
  TH2F *YZ_negativeX_hist_1= (TH2F*)my_file.Get("correction_dqdx_ZvsY_negativeX_hist_1");
  TH2F *YZ_positiveX_hist_1= (TH2F*)my_file.Get("correction_dqdx_ZvsY_positiveX_hist_1");
  ////////////////////// Importing Y-Z plane fractional corrections /////////////
  TH2F *YZ_negativeX_hist_0= (TH2F*)my_file.Get("correction_dqdx_ZvsY_negativeX_hist_0");
  TH2F *YZ_positiveX_hist_0= (TH2F*)my_file.Get("correction_dqdx_ZvsY_positiveX_hist_0");
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  TString output_file_str = "";
  if(useBirksModel) output_file_str = Form("./outputs/Xcalo_mich%s_r%d_%s.root",mn.Data(), run, ref_dedx_str.c_str());
  else output_file_str = Form("./outputs/Xcalo_mich%s_r%d_%s_Box.root",mn.Data(), run, ref_dedx_str.c_str());
  cout << "SB debug output_file_str : " << output_file_str << endl;
  //TFile *file = new TFile(Form("./outputs/Xcalo_mich%s_r%d_%s.root",mn.Data(), run, ref_dedx_str.c_str()),"recreate");
  TFile *file = new TFile(output_file_str, "recreate");
  TTree t1("t1","a simple Tree with simple variables");//creating a tree example
  Int_t run_number;
  Double_t event_time1;
  Float_t global_med_0,global_med_1,global_med_2; 
  t1.Branch("run_number",&run_number,"run_number/I");
  t1.Branch("event_time1",&event_time1,"event_time1/D");
  t1.Branch("global_med_0",&global_med_0,"global_med_0/F");
  t1.Branch("global_med_1",&global_med_1,"global_med_1/F");
  t1.Branch("global_med_2",&global_med_2,"global_med_2/F");
  Int_t runvalue=0;
  Double_t time1=0;
  //Filling the TTree

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t real_nentries = fChain->GetEntries();
  if (real_nentries >200000){
    cout<<"Total entries = "<<real_nentries<<endl;
    cout<<"Only use 200000 events."<<endl;
    nentries = 200000;
    real_nentries = nentries;
  }
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // == SB debug, original line
    //for (Long64_t jentry=0; jentry<1;jentry++) { // == SB debug
    Long64_t ientry = LoadTree(jentry);
    //cout << "SB debug, ref_dedx : " << ref_dedx << endl;

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%10000==0) cout<<jentry<<"/"<<real_nentries<<endl;
    if(jentry==0){
      time1=evttime;
      runvalue=run;
    }

    int x_bin;
    for(int i=0; i<cross_trks; ++i){
      bool testneg=0;
      bool testpos=0;
      if(trkstartx[i]*trkendx[i]>0) continue;
      if(trkstartx[i]<-350 or trkendx[i]<-350) testneg=1;      
      if(trkstartx[i]>350 or trkendx[i]>350) testpos=1;      
      //plane 2
      if(!((TMath::Abs(trkstartx[i])>350||trkstarty[i]<50||trkstarty[i]>550||trkstartz[i]<50||trkstartz[i]>645)&&(TMath::Abs(trkendx[i])>350||trkendy[i]<50||trkendy[i]>550||trkendz[i]<50||trkendz[i]>645))) continue;
      //  if(!(((TMath::Abs(trackthetaxz[i])>1.13) && (TMath::Abs(trackthetaxz[i])<2.0))||(TMath::Abs(trackthetayz[i])>1.22 && TMath::Abs(trackthetayz[i])<1.92))){
      if(!((abs(180/TMath::Pi()*trackthetaxz[i])>60 && abs(180/TMath::Pi()*trackthetaxz[i])<120)||(abs(180/TMath::Pi()*trackthetayz[i])>80 && abs(180/TMath::Pi()*trackthetayz[i])<100))){
	for(int j=1; j<TMath::Min(ntrkhits[i][2]-1,3000); ++j){
	  if((trkhity[i][2][j]<600)&&(trkhity[i][2][j]>0)){
	    if((trkhitz[i][2][j]<695)&&(trkhitz[i][2][j]>0)){
	      if(trkhitx[i][2][j]<0 && trkhitx[i][2][j]>-360 && testneg){//negative drift
		if(trkhitx[i][2][j]<0 && trkhitx[i][2][j+1]>0) continue;
		if(trkhitx[i][2][j]<0 && trkhitx[i][2][j-1]>0) continue;
		x_bin=dqdx_X_hist_2->FindBin(trkhitx[i][2][j]);
		float YZ_correction_factor_negativeX_2=YZ_negativeX_hist_2->GetBinContent(YZ_negativeX_hist_2->FindBin(trkhitz[i][2][j],trkhity[i][2][j]));
		float recom_correction = 1.;
		if(useBirksModel) {recom_correction = recom_factor_Birks(tot_Ef(trkhitx[i][2][j],trkhity[i][2][j],trkhitz[i][2][j]));}
		else{ recom_correction = recom_factor(tot_Ef(trkhitx[i][2][j],trkhity[i][2][j],trkhitz[i][2][j]));}
		float corrected_dqdx_2=trkdqdx[i][2][j]*YZ_correction_factor_negativeX_2*recom_correction;
		if (mn!="3") dqdx_value_2[x_bin-1].push_back(corrected_dqdx_2);
                else dqdx_value_2[x_bin-1].push_back(trkdqdx[i][2][j]);
		dedx_value_2[x_bin-1].push_back(trkdedx[i][2][j]);
                hdqdx[2]->Fill(trkdqdx[i][2][j]);
                hdedx[2]->Fill(trkdedx[i][2][j]);
	      }//X containment
	      if(trkhitx[i][2][j]>0 && trkhitx[i][2][j]<360 && testpos){//positive drift
		if(trkhitx[i][2][j]>0 && trkhitx[i][2][j+1]<0) continue;
		if(trkhitx[i][2][j]>0 && trkhitx[i][2][j-1]<0) continue;
		x_bin=dqdx_X_hist_2->FindBin(trkhitx[i][2][j]);
		float YZ_correction_factor_positiveX_2=YZ_positiveX_hist_2->GetBinContent(YZ_positiveX_hist_2->FindBin(trkhitz[i][2][j],trkhity[i][2][j]));
		float recom_correction = 1.;
		if(useBirksModel) {recom_correction = recom_factor_Birks(tot_Ef(trkhitx[i][2][j],trkhity[i][2][j],trkhitz[i][2][j]));}
                else{ recom_correction = recom_factor(tot_Ef(trkhitx[i][2][j],trkhity[i][2][j],trkhitz[i][2][j]));}
		float corrected_dqdx_2=trkdqdx[i][2][j]*YZ_correction_factor_positiveX_2*recom_correction;
                //std::cout<<TrkID[i]<<" "<<trkhitx[i][2][j]<<" "<<trkhity[i][2][j]<<" "<<trkhitz[i][2][j]<<" "<<YZ_correction_factor_positiveX_2<<" "<<trkdqdx[i][2][j]<<" "<<corrected_dqdx_2<<std::endl;
		if (mn!="3") dqdx_value_2[x_bin-1].push_back(corrected_dqdx_2);
                else dqdx_value_2[x_bin-1].push_back(trkdqdx[i][2][j]);
		dedx_value_2[x_bin-1].push_back(trkdedx[i][2][j]);
                hdqdx[2]->Fill(trkdqdx[i][2][j]);
                hdedx[2]->Fill(trkdedx[i][2][j]);
	      }//X containment
	    } // Z containment
	  } // Y containment
	} // loop over hits of the track in the given plane
      }//angular cut
    

      //// == Angular cut for plane 0 and 1
      //if(((abs(180/TMath::Pi()*trackthetaxz[i])>60 && abs(180/TMath::Pi()*trackthetaxz[i])<120))) continue;

      ////plane_1
      for(int j=1; j<TMath::Min(ntrkhits[i][1]-1,3000); ++j){
	if((trkhity[i][1][j]<600)&&(trkhity[i][1][j]>0)){
	  if((trkhitz[i][1][j]<695)&&(trkhitz[i][1][j]>0)){
	    if(trkhitx[i][1][j]<0 && trkhitx[i][1][j]>-360 && testneg){
	      if(abs(180/TMath::Pi()*trackthetaxz[i])>140){
	    	if(trkhitx[i][1][j]<0 && trkhitx[i][1][j+1]>0) continue;
		if(trkhitx[i][1][j]<0 && trkhitx[i][1][j-1]>0) continue;
		x_bin=dqdx_X_hist_1->FindBin(trkhitx[i][1][j]);
		float recom_correction=1.;
	     	if(useBirksModel) {recom_correction = recom_factor_Birks(tot_Ef(trkhitx[i][1][j],trkhity[i][1][j],trkhitz[i][1][j]));}
                else{ recom_correction = recom_factor(tot_Ef(trkhitx[i][1][j],trkhity[i][1][j],trkhitz[i][1][j]));}
		float YZ_correction_factor_negativeX_1=YZ_negativeX_hist_1->GetBinContent(YZ_negativeX_hist_1->FindBin(trkhitz[i][1][j],trkhity[i][1][j]));
		float corrected_dqdx_1=trkdqdx[i][1][j]*YZ_correction_factor_negativeX_1*recom_correction;
		if (mn!="3") dqdx_value_1[x_bin-1].push_back(corrected_dqdx_1);
                else dqdx_value_1[x_bin-1].push_back(trkdqdx[i][1][j]);
		dedx_value_1[x_bin-1].push_back(trkdedx[i][1][j]);
                hdqdx[1]->Fill(trkdqdx[i][1][j]);
                hdedx[1]->Fill(trkdedx[i][1][j]);
	      }
	    }
	    if(trkhitx[i][1][j]>0 && trkhitx[i][1][j]<360 && testpos){
	      if(trkhitx[i][1][j]>0 && trkhitx[i][1][j+1]<0) continue;
	      if(trkhitx[i][1][j]>0 && trkhitx[i][1][j-1]<0) continue;
	      if(abs(180/TMath::Pi()*trackthetaxz[i])<40){
		x_bin=dqdx_X_hist_1->FindBin(trkhitx[i][1][j]);
		float recom_correction=1.;
		if(useBirksModel) {recom_correction = recom_factor_Birks(tot_Ef(trkhitx[i][1][j],trkhity[i][1][j],trkhitz[i][1][j]));}
                else{ recom_correction = recom_factor(tot_Ef(trkhitx[i][1][j],trkhity[i][1][j],trkhitz[i][1][j]));}
		float YZ_correction_factor_positiveX_1=YZ_positiveX_hist_1->GetBinContent(YZ_positiveX_hist_1->FindBin(trkhitz[i][1][j],trkhity[i][1][j]));
		float corrected_dqdx_1=trkdqdx[i][1][j]*YZ_correction_factor_positiveX_1*recom_correction;
		if (mn!="3") dqdx_value_1[x_bin-1].push_back(corrected_dqdx_1);
                else dqdx_value_1[x_bin-1].push_back(trkdqdx[i][1][j]);
		dedx_value_1[x_bin-1].push_back(trkdedx[i][1][j]);
                hdqdx[1]->Fill(trkdqdx[i][1][j]);
                hdedx[1]->Fill(trkdedx[i][1][j]);
	      }
	    }
	  } // Z containment
	} // Y containment
      } // loop over hits of the track in the given plane
   
      /////plane_0
      for(int j=1; j<TMath::Min(ntrkhits[i][0]-1,3000); ++j){
	if((trkhity[i][0][j]<600)&&(trkhity[i][0][j]>0)){
	  if((trkhitz[i][0][j]<695)&&(trkhitz[i][0][j]>0)){
	    if(trkhitx[i][0][j]<0 && trkhitx[i][0][j]>-360 && testneg){
	      if(abs(180/TMath::Pi()*trackthetaxz[i])<40){
		if(trkhitx[i][0][j]<0 && trkhitx[i][0][j+1]>0) continue;
		if(trkhitx[i][0][j]<0 && trkhitx[i][0][j-1]>0) continue;
		x_bin=dqdx_X_hist_0->FindBin(trkhitx[i][0][j]);
		float recom_correction=1.;
		if(useBirksModel) {recom_correction = recom_factor_Birks(tot_Ef(trkhitx[i][0][j],trkhity[i][0][j],trkhitz[i][0][j]));}
                else{ recom_correction = recom_factor(tot_Ef(trkhitx[i][0][j],trkhity[i][0][j],trkhitz[i][0][j]));}
		float YZ_correction_factor_negativeX_0=YZ_negativeX_hist_0->GetBinContent(YZ_negativeX_hist_0->FindBin(trkhitz[i][0][j],trkhity[i][0][j]));
		float corrected_dqdx_0=trkdqdx[i][0][j]*YZ_correction_factor_negativeX_0*recom_correction;
		if (mn!="3") dqdx_value_0[x_bin-1].push_back(corrected_dqdx_0);
                else  dqdx_value_0[x_bin-1].push_back(trkdqdx[i][0][j]);
		dedx_value_0[x_bin-1].push_back(trkdedx[i][0][j]);
                hdqdx[0]->Fill(trkdqdx[i][0][j]);
                hdedx[0]->Fill(trkdedx[i][0][j]);
                //if (x_bin == 1) cout<<event<<" neg "<<x_bin<<" "<<trkhitx[i][0][j]<<" "<<trkhity[i][0][j]<<" "<<trkhitz[i][0][j]<<endl;
	      }
	    }
	    if(trkhitx[i][0][j]>0 && trkhitx[i][0][j]<360 && testpos){
	      if(abs(180/TMath::Pi()*trackthetaxz[i])>140){
		if(trkhitx[i][0][j]>0 && trkhitx[i][0][j+1]<0) continue;
		if(trkhitx[i][0][j]>0 && trkhitx[i][0][j-1]<0) continue;
		x_bin=dqdx_X_hist_0->FindBin(trkhitx[i][0][j]);
		float recom_correction=1.;
		if(useBirksModel) {recom_correction = recom_factor_Birks(tot_Ef(trkhitx[i][0][j],trkhity[i][0][j],trkhitz[i][0][j]));}
                else{ recom_correction = recom_factor(tot_Ef(trkhitx[i][0][j],trkhity[i][0][j],trkhitz[i][0][j]));}
		float YZ_correction_factor_positiveX_0=YZ_positiveX_hist_0->GetBinContent(YZ_positiveX_hist_0->FindBin(trkhitz[i][0][j],trkhity[i][0][j]));
		float corrected_dqdx_0=trkdqdx[i][0][j]*YZ_correction_factor_positiveX_0*recom_correction;
		if (mn!="3") dqdx_value_0[x_bin-1].push_back(corrected_dqdx_0);
                else  dqdx_value_0[x_bin-1].push_back(trkdqdx[i][0][j]);
		dedx_value_0[x_bin-1].push_back(trkdedx[i][0][j]);
                hdqdx[0]->Fill(trkdqdx[i][0][j]);
                hdedx[0]->Fill(trkdedx[i][0][j]);
                //if (x_bin == 1) cout<<event<<" pos "<<x_bin<<" "<<trkhitx[i][0][j]<<" "<<trkhity[i][0][j]<<" "<<trkhitz[i][0][j]<<endl;
	      }
	    }
	  } // Z containment
	} // Y containment
      } // loop over hits of the track in the given plane

    } // loop over crossing tracks in the event
  } // loop over jentries

  std::cout << "*************** Calculating the local median dQ/dx values for each Y-Z cell ******************" << std::endl;

  ////////////////////// Getting inforamtion from dqdx_value vector ////////////////////
  for(size_t i=0; i<dqdx_value_2.size(); i++){
    //std::cout<<"no of entries in each bin plane 2 "<<i<<"  "<<dqdx_value_2[i].size()<<std::endl;
    if(dqdx_value_2[i].size()>5){
      for(size_t k=0; k<dqdx_value_2[i].size(); k++){
	all_dqdx_value_2.push_back(dqdx_value_2[i][k]);
      }
      float local_median_dqdx_2=TMath::Median(dqdx_value_2[i].size(),&dqdx_value_2[i][0]);
      dqdx_X_hist_2->SetBinContent(i+1,local_median_dqdx_2);
      float local_median_dedx_2=TMath::Median(dedx_value_2[i].size(),&dedx_value_2[i][0]);
      dedx_X_hist_2->SetBinContent(i+1,local_median_dedx_2);
    }
  }
  if (run>10000) run = 0;
  float global_median_dqdx_2=TMath::Median(all_dqdx_value_2.size(),&all_dqdx_value_2[0]); 
  global_med_2=global_median_dqdx_2;//Filling the Tree variable
  ofstream outfile0,outfile1,outfile2;
  outfile2.open(Form("./outputs/global_median_2_r%d_%s.txt", run, ref_dedx_str.c_str()));
  outfile2<<run<<"\t"<<global_median_dqdx_2<<std::endl;
  //////////////////////////////////////////////////////////////////////////////////////
 
  std::cout << "**************** Calculating fractional correction for each x cell *********************" << std::endl;
 
  ///////////////////////// Calculating fractional corrections /////////////////////////////
 
  for(size_t i=0; i<dqdx_value_2.size(); i++){
    if(dqdx_value_2[i].size()>5){
      float local_median_dqdx_2=TMath::Median(dqdx_value_2[i].size(),&dqdx_value_2[i][0]);
      float fractional_dqdx_2=float(global_median_dqdx_2)/local_median_dqdx_2;
      dqdx_X_correction_hist_2->SetBinContent(i+1,fractional_dqdx_2);
      dqdx_frac_correction_2[i].push_back(fractional_dqdx_2);
    }
  }
 
  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << "**************** Calculating XYZ corrected dQ/dx values ********************" << std::endl;
  //////////////// How corrected X dqdx distribution looks like /////////////////////
  for(size_t i=0; i<dqdx_value_2.size(); i++){
    if(dqdx_value_2[i].size()>5){
      float local_median_dqdx_2=TMath::Median(dqdx_value_2[i].size(),&dqdx_value_2[i][0]);
      float Xcorrected_dqdx_2=local_median_dqdx_2*dqdx_frac_correction_2[i][0];
      corrected_dqdx_X_hist_2->SetBinContent(i+1,Xcorrected_dqdx_2);
    }
  }
 
  //////////////////////////////////////////////////////////////////////////////////
 
  dqdx_X_hist_2->Write();
  dedx_X_hist_2->Write();
  dqdx_X_correction_hist_2->Write();
  corrected_dqdx_X_hist_2->Write();
 ////////////////////// Getting inforamtion from dqdx_value vector ////////////////////
  for(size_t i=0; i<dqdx_value_1.size(); i++){
    //std::cout<<"no of entries in each bin plane 1 "<<i<<"  "<<dqdx_value_1[i].size()<<std::endl;
    if(dqdx_value_1[i].size()>5){
      for(size_t k=0; k<dqdx_value_1[i].size(); k++){
	all_dqdx_value_1.push_back(dqdx_value_1[i][k]);
      }
      float local_median_dqdx_1=TMath::Median(dqdx_value_1[i].size(),&dqdx_value_1[i][0]);
      dqdx_X_hist_1->SetBinContent(i+1,local_median_dqdx_1);
      float local_median_dedx_1=TMath::Median(dedx_value_1[i].size(),&dedx_value_1[i][0]);
      dedx_X_hist_1->SetBinContent(i+1,local_median_dedx_1);
    }
  }
 
  float global_median_dqdx_1=TMath::Median(all_dqdx_value_1.size(),&all_dqdx_value_1[0]);
  global_med_1=global_median_dqdx_1;//Filling the Tree variable
  outfile1.open(Form("./outputs/global_median_1_r%d_%s.txt", run, ref_dedx_str.c_str()));
  outfile1<<run<<"\t"<<global_median_dqdx_1<<std::endl; 
 
  //////////////////////////////////////////////////////////////////////////////////////
 
  std::cout << "**************** Calculating fractional correction for each x cell *********************" << std::endl;
 
  ///////////////////////// Calculating fractional corrections /////////////////////////////
 
  for(size_t i=0; i<dqdx_value_1.size(); i++){
    if(dqdx_value_1[i].size()>5){
      float local_median_dqdx_1=TMath::Median(dqdx_value_1[i].size(),&dqdx_value_1[i][0]);
      float fractional_dqdx_1=float(global_median_dqdx_1)/local_median_dqdx_1;
      dqdx_X_correction_hist_1->SetBinContent(i+1,fractional_dqdx_1);
      dqdx_frac_correction_1[i].push_back(fractional_dqdx_1);
    }
  }
 
  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << "**************** Calculating XYZ corrected dQ/dx values ********************" << std::endl;
  //////////////// How corrected X dqdx distribution looks like /////////////////////
  for(size_t i=0; i<dqdx_value_1.size(); i++){
    if(dqdx_value_1[i].size()>5){
      float local_median_dqdx_1=TMath::Median(dqdx_value_1[i].size(),&dqdx_value_1[i][0]);
      float Xcorrected_dqdx_1=local_median_dqdx_1*dqdx_frac_correction_1[i][0];
      corrected_dqdx_X_hist_1->SetBinContent(i+1,Xcorrected_dqdx_1);
    }
  }
  //////////////////////////////////////////////////////////////////////////////////
  dqdx_X_hist_1->Write();
  dedx_X_hist_1->Write();
  dqdx_X_correction_hist_1->Write();
  corrected_dqdx_X_hist_1->Write();


 ////////////////////// Getting inforamtion from dqdx_value vector ////////////////////
  for(size_t i=0; i<dqdx_value_0.size(); i++){
    //std::cout<<"no of entries in each bin plane 0 "<<i<<"  "<<dqdx_value_0[i].size()<<std::endl;
    if(dqdx_value_0[i].size()>5){
      for(size_t k=0; k<dqdx_value_0[i].size(); k++){
	all_dqdx_value_0.push_back(dqdx_value_0[i][k]);
      }
      float local_median_dqdx_0=TMath::Median(dqdx_value_0[i].size(),&dqdx_value_0[i][0]);
      dqdx_X_hist_0->SetBinContent(i+1,local_median_dqdx_0);
      float local_median_dedx_0=TMath::Median(dedx_value_0[i].size(),&dedx_value_0[i][0]);
      dedx_X_hist_0->SetBinContent(i+1,local_median_dedx_0);
    }
  }
  float global_median_dqdx_0=TMath::Median(all_dqdx_value_0.size(),&all_dqdx_value_0[0]); 
  global_med_0=global_median_dqdx_0;//Filling the Tree variable
  outfile0.open(Form("./outputs/global_median_0_r%d_%s.txt",run, ref_dedx_str.c_str()));
  outfile0<<run<<"\t"<<global_median_dqdx_0<<std::endl; 
  run_number=runvalue;
  event_time1=time1;
  t1.Fill();//Filling the Tree
  //////////////////////////////////////////////////////////////////////////////////////
 
  std::cout << "**************** Calculating fractional correction for each x cell *********************" << std::endl;
 
  ///////////////////////// Calculating fractional corrections /////////////////////////////
 
  for(size_t i=0; i<dqdx_value_0.size(); i++){
    if(dqdx_value_0[i].size()>5){
      float local_median_dqdx_0=TMath::Median(dqdx_value_0[i].size(),&dqdx_value_0[i][0]);
      float fractional_dqdx_0=float(global_median_dqdx_0)/local_median_dqdx_0;
      dqdx_X_correction_hist_0->SetBinContent(i+1,fractional_dqdx_0);
      dqdx_frac_correction_0[i].push_back(fractional_dqdx_0);
    }
  }
 
  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << "**************** Calculating XYZ corrected dQ/dx values ********************" << std::endl;
  //////////////// How corrected X dqdx distribution looks like /////////////////////
  for(size_t i=0; i<dqdx_value_0.size(); i++){
    if(dqdx_value_0[i].size()>5){
      float local_median_dqdx_0=TMath::Median(dqdx_value_0[i].size(),&dqdx_value_0[i][0]);
      float Xcorrected_dqdx_0=local_median_dqdx_0*dqdx_frac_correction_0[i][0];
      corrected_dqdx_X_hist_0->SetBinContent(i+1,Xcorrected_dqdx_0);
    }
  }
 
  //////////////////////////////////////////////////////////////////////////////////
 
  dqdx_X_hist_0->Write();
  dedx_X_hist_0->Write();
  dqdx_X_correction_hist_0->Write();
  corrected_dqdx_X_hist_0->Write();

  for (int i = 0; i<3; ++i){
    hdqdx[i]->Write();
    hdedx[i]->Write();
  }

  file->Close(); 
  dqdx_X_hist_2->Draw();
  TFile treefile(Form("./outputs/globalmedians_cathanode_r%d_%s.root",run, ref_dedx_str.c_str()),"RECREATE");
  t1.Write();


  std::cout << "*************** X_Correction_make_class.C macro has ended ******************" << std::endl; 
}

int main(int argc, char *argv[]) {
  
  if (!argv[2]) {
    cout << "Error: No input file or michelremoving tree number was provided!" << endl;
    cout << "Usage: " << endl;
    cout << "make_x_correction root_file_or_list [michelremoving_tree_number]" << endl;

    return 0;
  }
  
  string infile = argv[1];
  string michelnumber = argv[2];
  string sce = argv[3];
  ref_dedx_str = argv[4];
  
  if (!(michelnumber == "0"||michelnumber == "1"||michelnumber == "2"||michelnumber == "3")){
    cout << "Error: Michel tree number must be 0,1, or 2" << endl;
    return 0;
    }

  if (michelnumber=="0") michelnumber = "";
  cout << Form("michelremoving%s/Event", michelnumber.c_str()) << endl;

  TChain* shtree = new TChain("Event");

  if (infile.substr(infile.find_last_of(".") + 1) == "root"){
    shtree->Add(Form("%s/michelremoving%s/Event", infile.c_str(), michelnumber.c_str()));
  }

  else /*if(infile.substr(infile.find_last_of(".") + 1) == "list" || infile.substr(infile.find_last_of(".") + 1) == "txt")*/{
    std::ifstream in;
    in.open(infile.c_str());
    char line[1024];

    while(1){
      in.getline(line,1024);
      if (!in.good()) break;
      shtree->Add(Form("%s/michelremoving%s/Event", line, michelnumber.c_str()));
    }
    in.close();
    in.clear();
  }

  if (sce=="0"){
    cout<<"SCE off"<<endl;
    sceon = false;
    userecom = false;
  }
  else{
    sceon = true;
    userecom = true;
    cout<<"SCE on"<<endl;
  }

  ref_dedx = stof(ref_dedx_str);
  cout << "Reference dE/dx : " << ref_dedx << endl;

  protoDUNE_X_calib_Birks t(shtree);
  
  t.Loop(michelnumber.c_str());
} // main
