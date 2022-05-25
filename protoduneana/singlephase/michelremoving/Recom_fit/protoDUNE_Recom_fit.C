#define protoDUNE_Recom_fit_cxx
#include "protoduneana/singlephase/michelremoving/scripts/LanGausFit.h"
#include "protoDUNE_Recom_fit.h"
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
#include <TSpline.h>
#include <TText.h>
#include <TFrame.h>
#include <TMinuit.h>
#include <TVectorD.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace std;

//int hitplane=1;
//bool usemap=true;
//const double calib_factor =4.50e-3; // ******************* change the calibraion factor here with the normalisation ****************** //
//double calib_factor =1.029e-3;
double pitchvalue=0.65;
//double Efield = 0.50;//kV/cm protoDUNE electric filed
//double normalisation_factor[3]={0.992,0.990,0.989};//for plane 0
const int Z=18; //Atomic number of Argon
const double A=39.948; // g/mol Atomic mass of Argon
const double I=188.0e-6; // ev
const double K=0.307; // Mev.cm^2 / mol
const double Mmu=105.658; // Mev for Mu
const double Me=0.51; // Mev for electron
const double rho=1.396;//g/cm^3
string outfile_name;
string sce;
TString mn = "2";
bool use_Birk_model = true;

double spline_KE[13];
double spline_Range[13];
//TFile * OpenFile(const std::string filename);
TFile * OpenFile(const std::string filename) {
  TFile * theFile = 0x0;
  mf::LogInfo("OpenFile") << "Searching for " << filename;
  if (cet::file_exists(filename)) {
    mf::LogInfo("OpenFile") << "File exists. Opening " << filename;
    theFile = new TFile(filename.c_str());
    if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
      delete theFile;
      theFile = 0x0;
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not open " << filename;
    }
  }
  else {
    mf::LogInfo("OpenFile") << "File does not exist here. Searching FW_SEARCH_PATH";
    cet::search_path sp{"FW_SEARCH_PATH"};
    std::string found_filename;
    auto found = sp.find_file(filename, found_filename);
    if (!found) {
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not find " << filename;
    }

    mf::LogInfo("OpenFile") << "Found file " << found_filename;
    theFile = new TFile(found_filename.c_str());
    if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
      delete theFile;
      theFile = 0x0;
      throw cet::exception("ProtoDUNECalibration.cxx") << "Could not open " << found_filename;
    }
  }
  return theFile;
}

////getting the variable Efield using data driven maps
//TFile *ef/*=new TFile("SCE_DataDriven_180kV_v4.root")*/;
//TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
//TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
//TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
//TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
//TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
//TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
float protoDUNE_Recom_fit::tot_Ef(float xval,float yval,float zval){
  float E0value=0.4867;
  //if(!usemap) return E0value;
  //std::cout << xval << " " << yval << " " << zval << std::endl;
  if (sce=="on"){
    if(xval>=0){
      for (auto h : pos_hists) {
        if (h->GetXaxis()->FindBin(xval) < 1 ||
            h->GetXaxis()->FindBin(xval) > h->GetNbinsX()) {
          std::cout << "xval oob: " << xval << std::endl;
        }
        if (h->GetYaxis()->FindBin(yval) < 1 ||
            h->GetYaxis()->FindBin(yval) > h->GetNbinsY()) {
          std::cout << "yval oob: " << yval << std::endl;
        }
        if (h->GetZaxis()->FindBin(zval) < 1 ||
            h->GetZaxis()->FindBin(zval) > h->GetNbinsZ()) {
          std::cout << "zval oob: " << zval << std::endl;
        }
      }
      float ex=E0value+E0value*xpos->Interpolate(xval,yval,zval);
      float ey=E0value*ypos->Interpolate(xval,yval,zval);
      float ez=E0value*zpos->Interpolate(xval,yval,zval);
      return sqrt(ex*ex+ey*ey+ez*ez);
    }
    if(xval<0){
      for (auto h : neg_hists) {
        if (h->GetXaxis()->FindBin(xval) < 1 ||
            h->GetXaxis()->FindBin(xval) > h->GetNbinsX()) {
          std::cout << "xval oob: " << xval << std::endl;
        }
        if (h->GetYaxis()->FindBin(yval) < 1 ||
            h->GetYaxis()->FindBin(yval) > h->GetNbinsY()) {
          std::cout << "yval oob: " << yval << std::endl;
        }
        if (h->GetZaxis()->FindBin(zval) < 1 ||
            h->GetZaxis()->FindBin(zval) > h->GetNbinsZ()) {
          std::cout << "zval oob: " << zval << std::endl;
        }
      }
      float ex=E0value+E0value*xneg->Interpolate(xval,yval,zval);
      float ey=E0value*yneg->Interpolate(xval,yval,zval);
      float ez=E0value*zneg->Interpolate(xval,yval,zval);
      return sqrt(ex*ex+ey*ey+ez*ez);
    }
  }
  return E0value;
}

////********************************************///


double beta(double gamma){
  double value=TMath::Sqrt(1-(1.0/(gamma*gamma)));
  return value;
}

double gamma(double KE,double mass){
  double value=(double(KE)/mass)+1;
  return value;
}

double KE=266.;
double g=gamma(KE,Mmu);
double b=beta(g);

const double C=-5.2146;
const double X0=0.2;
const double X1=3.0;
const double a=0.19559;
const double m=3.0;
const double N=2*TMath::Log(10);

double density(double bg){//replaced x by x1
  double value;
  double x = TMath::Log10(bg);
  if(x<X0) return 0;
  if(x>X1) return N*x + C;
  value=a*(TMath::Power((X1-x),m));
  return N*x + C + value;
}

double Wmax(double KE,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double num=2*Me*(TMath::Power(b*g,2));
  double den=1+double(2*g*Me)/mass + TMath::Power((double(Me)/mass),2);
  double value=double(num)/den;
  return value;
}

double dEdx(double KE,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double f=K*(double(Z)/A)*(TMath::Power(1.0/b,2));
  double wmax=Wmax(KE,mass);
  double a0=0.5*(TMath::Log(double(2*Me*(TMath::Power(b*g,2))*wmax)/(I*I)));
  double dens=density(b*g);
  double value=f*rho*(a0-b*b-double(dens)/2);
  return value;
}

double dpdx(double KE,double x,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double epsilon=(double(K)/2)*(double(Z)/A)*(double(x*rho)/(b*b));
  double A0=double(2*Me*(TMath::Power((b*g),2)))/I;
  double A1=double(epsilon)/I;
  double value=(1.0/x)*epsilon*((TMath::Log(A0)) + TMath::Log(A1) + 0.2 - b*b - density(b*g));
  return value;
}

///////////////////// End of Landau-Vavilov function /////////////////////////////////////////////

///////////////////////////////////// Function definition ////////////////////////////////

Double_t langaufun(Double_t *x, Double_t *par) {
  Double_t invsq2pi = 0.398942280401;// Control constants
  //Double_t mpshift = -0.22278298;
  Double_t np = 500.0;
  Double_t sc = 5.0;// convolution extends to +-sc Gaussian sigmas 
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  //mpc = par[1]- mpshift * par[0];
  mpc=par[1];
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow)/np;
 
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

//////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////// Function definition ////////////////////////////////

 TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status)
// TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  Int_t i;
  Char_t FunName[100];

  sprintf(FunName,"Fitfcn_%s",his->GetName());

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MPV","Area","GSigma");
   
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

//  his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
  TFitResultPtr fitres = his->Fit(FunName,"RBOSQ"); // fit within specified range, use ParLimits, do not plot /////////////////// Initial code use the mode "RBO" (commented by VARUNA) ///////////

  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
 
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  Status[0] = fitres->CovMatrixStatus();
 
  return (ffit);              // return fit function
}

/////////////////////////////// Function definition /////////////////////////////////////

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
  Double_t p,x,fy,fxr,fxl;
  Double_t step;
  Double_t l,lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;

  // Search for maximum

  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l = -1.0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = langaufun(&x,params);
    if (l < lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-1);

  maxx = x;
  fy = l/2;

  // Search for right x location of fy

  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-2);

  fxr = x;

  // Search for left x location of fy

  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-3);

  fxl = x;
  FWHM = fxr - fxl;
  return (0);
}

////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////// Function definition ////////////////////////////////

double Wion = 23.6e-6; //parameter from ArgoNeuT experiment at 0.481kV/cm
Float_t Dedx(float dqdx, float Ef){
  double Rho = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
  //double Wion = 23.6e-6;//parameter from ArgoNeuT experiment at 0.481kV/cm
  //ArgoNeuT parameters
  double alpha = 0.93;//parameter from ArgoNeuT experiment at 0.481kV/cm 
  double betap = 0.212;//(kV/cm)(g/cm^2)/MeV
  //Abbey's parameters https://indico.fnal.gov/event/46502/contributions/206721/attachments/139268/174710/recombination_20210126.pdf
  //double alpha = 0.854;//parameter from ArgoNeuT experiment at 0.481kV/cm 
  //double betap = 0.208;//(kV/cm)(g/cm^2)/MeV
  //Abbey's parameters https://indico.fnal.gov/event/46503/contributions/215375/attachments/143230/181102/recombination_20210519.pdf
//  double alpha = 0.912;//parameter from ArgoNeuT experiment at 0.481kV/cm 
//  double betap = 0.195;//(kV/cm)(g/cm^2)/MeV

  //cout << "SB debug [Dedx] dqdx : " << dqdx << ", Ef : " << Ef << endl;

  return (exp(dqdx*(betap/(Rho*Ef)*Wion))-alpha)/(betap/(Rho*Ef));
}

Float_t Dedx_Birk(float dqdx, float Ef){
  double Rho = 1.383; // == g/cm^3 (liquid argon density at a pressure 18.0 psia)
  //double Wion = 23.6e-6; // == MeV from ArgoNeuT experiment at 0.481kV/cm
  //double A_B = 0.806; //  == ArgoNeuT :  at 0.481kV/cm
  //double k_B = 0.052; // == ArgoNeuT (kV/cm)(g/cm2 )/MeV
  double A_B = 0.800; // == ICARUS : at 200, 350, 500 V/cm
  double k_B = 0.0486; // == ICARUS (kV/cm)(g/cm2 )/MeV 

  return( dqdx / ( (A_B/Wion) - k_B * (dqdx /(Rho * Ef) )) );

}

//////////////////////////////////////////////////////////////////////////////////////

void protoDUNE_Recom_fit::LoopLite(std::vector<double> & norm_factors,
                                    vector<double> & CalibConstants,
                                    TFile & outfile) {

  // == vector<vector<double>> & calib_factors,
  std::cout << "******************************* LoopLite is running *******************************" << std::endl;

  outfile.cd();

  size_t nbin=12; // == SB debug, (dE/dx_theory) / Etot variable number of bins
  double binsize=0.05; // == SB debug, (dE/dx_theory) / Etot variable bin size
  double bin_start = 3.0; // == SB debug, (dE/dx_theory) / Etot variable starts from
  //TH1D *dedx[nbin];
  //{hitplane, {calib_index, {detector bins}}}
  std::map<int, std::vector<TH1D*>> dedx_measured_hist;
  std::map<int, std::vector<TH1D*>> dedx_theory_hist;

  vector<vector<vector<float>>> dedx_measured_value(3);
  vector<vector<vector<float>>> dedx_theory_value(3);

  vector<vector<float>> dedx_measured_median(3);
  vector<vector<float>> dedx_theory_median(3);
  vector<vector<float>> dedx_measured_median_err(3);
  vector<vector<float>> dedx_theory_median_err(3);

  vector<vector<float>> Recom_factor(3);
  vector<vector<float>> Recom_factor_err(3);

  for(size_t i = 0; i < 3; i++) {
    dedx_measured_value[i].resize(12);
    dedx_theory_value[i].resize(12);

    dedx_measured_median[i].resize(12);
    dedx_theory_median[i].resize(12);
    dedx_measured_median_err[i].resize(12);
    dedx_theory_median_err[i].resize(12);
    Recom_factor[i].resize(12);
    Recom_factor_err[i].resize(12);

    //dedx_measured_hist[i].push_back(std::vector<TH1D*>());
    //dedx_theory_hist[i].push_back(std::vector<TH1D*>());
    dedx_measured_hist[i] = std::vector<TH1D*>();
    dedx_theory_hist[i] = std::vector<TH1D*>();
    for(size_t j = 0; j < nbin; j++){
      dedx_measured_hist[i].push_back(new TH1D(Form("dedx_measured_%zu_%zu", i, j), Form("Plane:%zu %zu bin of dedx_theory/Etot", i, j), 200, 0., 10.));
      dedx_theory_hist[i].push_back(new TH1D(Form("dedx_theory_%zu_%zu", i, j), Form("Plane:%zu %zu bin of dedx_theory/Etot", i, j), 200, 0., 10.));
    }
  }

  //for (int i=0; i<nbin; ++i){
  //  if(i==0) dedx[i] = new TH1D(Form("dedx_%d",i),Form("dedx_%d",i),300,0.0,15); 
  //  if(i!=0) dedx[i] = new TH1D(Form("dedx_%d",i),Form("dedx_%d",i),200,0.0,10);
  //  dedx[i]->SetLineColor(kBlack); 
  //  dedx[i]->Sumw2();
  //}
 
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
 
  ///////////////// Make any changes to the Y and Z bin sizes here ///////////////

  //int x_bin_size = 5; // 148 bins in x direction
  fChain->GetEntry(0); 

  /////////////////////Importing X fractional corrections//////////////////////

  std::vector<TH1F*> X_correction_hists;
  std::vector<TH2F*> YZ_correction_neg_hists, YZ_correction_pos_hists;
  for (int i = 0; i < 3; ++i) {
    //TH1F *X_correction_hist = (TH1F*)fXFile->Get(Form("dqdx_X_correction_hist_%d",hitplane));
    //TH2F *YZ_correction_neg_hist=(TH2F*)fYZFile->Get(Form("correction_dqdx_ZvsY_negativeX_hist_%d",hitplane));
    //TH2F *YZ_correction_pos_hist=(TH2F*)fYZFile->Get(Form("correction_dqdx_ZvsY_positiveX_hist_%d",hitplane));
    X_correction_hists.push_back((TH1F*)fXFile->Get(Form("dqdx_X_correction_hist_%d",i)));
    YZ_correction_neg_hists.push_back((TH2F*)fYZFile->Get(Form("correction_dqdx_ZvsY_negativeX_hist_%d",i)));
    YZ_correction_pos_hists.push_back((TH2F*)fYZFile->Get(Form("correction_dqdx_ZvsY_positiveX_hist_%d",i)));
  }

  // == Histograms for theory dE/dx for all anode planes
  map<int, TH1D*> dedx_theory;
  map<int, TH1D*> dedx_measured;
  for(size_t i = 0; i< 3; i++){
    dedx_theory[i] = new TH1D(Form("dedx_theory_plane%zu", i), Form("Plane:%zu dEdx theory", i), 200, 0., 10.);
    dedx_measured[i] = new TH1D(Form("dedx_measured_plane%zu", i), Form("Plane:%zu dEdx measured", i), 60, 0., 3.);

  }

  TSpline3 *sp = new TSpline3("Cubic Spline", &spline_Range[0], &spline_KE[0],13,"b2e2",0,0);

  ////////////////////////////////////////////////////////////////////////////////// 

 
  if (fChain == 0) return;

  //double avgke[3] = {0};
  //double avgdx[3] = {0};
  //int nhits[3] = {0};

  Long64_t nentries = fChain->GetEntries();
  if (nentries > 200000){
    cout<<"Total entries = "<<nentries<<endl;
    cout<<"Only use 200000 events."<<endl;
    nentries = 200000;
  }
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<10000;jentry++) { // == SB debug

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if(!(jentry%10000)) cout<<jentry<<"/"<<nentries<< endl;
    for(int i=0; i<cross_trks; ++i){
      if (trkstartx[i]*trkendx[i]>0) continue;
      //if(peakT_min[i]<100 || peakT_max[i]>5900 ||
      //   trklen[i]<100 || trklen[i]>700 || (trkendz[i]>226 && trkendz[i]<236) ||
      //   (trkstartz[i]>226 && trkstartz[i]<236) ||
      //   (trkendz[i]>456 && trkendz[i]<472)||
      //   (trkstartz[i]>456 && trkstartz[i]<472)) continue;//filter for plane 2
      if(peakT_min[i]<100||peakT_max[i]>5900||trklen[i]<100||trklen[i]>700||(trkendz[i]>226 && trkendz[i]<236)||(trkstartz[i]>226 && trkstartz[i]<236)||(trkendz[i]>456 && trkendz[i]<472)||(trkstartz[i]>456 && trkstartz[i]<472)) continue;//filter for plane 2
      if(adjacent_hits[i]!=0 || dist_min[i]>5) continue;
      for (size_t ihitplane = 0; ihitplane < 3; ++ihitplane) {
        vector<float> res, dq, first5dq, last5dq;
        if(ihitplane==2 && ((abs(180/TMath::Pi()*trackthetaxz[i])>60 &&
                             abs(180/TMath::Pi()*trackthetaxz[i])<120)||
                            (abs(180/TMath::Pi()*trackthetayz[i])>80 &&
                             abs(180/TMath::Pi()*trackthetayz[i])<100))) continue;
 
        if(ihitplane == 2 && (lastwire[i]<=5 || lastwire[i]>=475)) continue;  //for plane 2  

        //TY: not sure if induction plane wires boundary is [0,800]
        if(ihitplane != 2 && (lastwire[i]<=5 || lastwire[i]>=795)) continue; 

        for(int j=0;j<ntrkhits[i][ihitplane];j++){
          res.push_back(trkresrange[i][ihitplane][j]);
          dq.push_back(trkdqdx[i][ihitplane][j]);
        }
        /**********************end of buffer filling*****************************/

        if(res.size()==0) continue;
        /***********************removed empty tracks to avoid segmentation fault******************/

        int siz1=ntrkhits[i][ihitplane];
        float max=*max_element(res.begin(),res.end());
       
        /************************flipping wrongly ordered residual range values****************************/
        bool test=true;
        if((trkhity[i][ihitplane][siz1-1]<trkhity[i][ihitplane][0] && trkresrange[i][ihitplane][siz1-1]>trkresrange[i][ihitplane][0])||(trkhity[i][ihitplane][0]<trkhity[i][ihitplane][siz1-1] && trkresrange[i][ihitplane][0]>trkresrange[i][ihitplane][siz1-1])){
          test=false;
          for(int i1=0;i1<ntrkhits[i][ihitplane];i1++){
            trkresrange[i][ihitplane][i1]=res[siz1-i1-1];
          }
          cout<<"This is a flipped track"<<endl;
        }
       
        /***************calculating the ratio of dQdx for first 5cm and last 5 cm of a track********************/ 
        for(size_t k=0;k<res.size();k++){
          if(trkresrange[i][ihitplane][k]<5) first5dq.push_back(dq[k]);
          if(trkresrange[i][ihitplane][k]>max-5) last5dq.push_back(dq[k]);
        }
       
        if(first5dq.size()<5){
          continue;
        }

        float med1 = TMath::Median(first5dq.size(), &first5dq[0]);
        float med2 = TMath::Median(last5dq.size(), &last5dq[0]);
        if(!((med1/med2)>1.4)) continue;
        if(!test) continue;

        for(int j=0; j<TMath::Min(ntrkhits[i][ihitplane],3000); ++j){
          float ke = sp->Eval(trkresrange[i][ihitplane][j]);
	  //cout << "SB debug trkresrange : " << trkresrange[i][ihitplane][j] << ", ke : " << ke << endl;
          if (ke<100 || ke>450) continue;
          if (trkpitch[i][ihitplane][j]>=0.5 && trkpitch[i][ihitplane][j]<=0.8 &&
              trkhity[i][ihitplane][j]>0 && trkhity[i][ihitplane][j]<600 &&
              trkhitz[i][ihitplane][j]>0 && trkhitz[i][ihitplane][j]<695) {

            TH2F * YZ_hist = 0x0;           

            if(trkhitx[i][ihitplane][j]>-360 && trkhitx[i][ihitplane][j]<0){ //negative X direction
              bool is_good = false;
              if(ihitplane == 1 && abs(180/TMath::Pi()*trackthetaxz[i])>140) {
                is_good = true;
              } //plane 1
              else if(ihitplane == 0 && abs(180/TMath::Pi()*trackthetaxz[i])<40){ //plane 0
                //std::cout << "skipping 1" << std::endl;
                is_good = true;
              }
              else if (ihitplane == 2) {
                is_good = true;
              }
              if (!is_good) continue;
              YZ_hist = YZ_correction_neg_hists[ihitplane];
            }
            else if(trkhitx[i][ihitplane][j]>0 && trkhitx[i][ihitplane][j]<360){ //positive X direction
              bool is_good = false;
              if(ihitplane == 1 && abs(180/TMath::Pi()*trackthetaxz[i])<40){ //plane 1
                is_good = true;
                //continue;
              }
              else if(ihitplane == 0 && abs(180/TMath::Pi()*trackthetaxz[i])>140){ //plane 0
                is_good = true;
                //std:!:cout << "skipping 3" << std::endl;
                //continue;
              }
              else if (ihitplane == 2) {
                is_good = true;
              }
              if (!is_good) continue;

              YZ_hist = YZ_correction_pos_hists[ihitplane];
            }
            else {
              continue;
            }

            int x_bin = X_correction_hists[ihitplane]->FindBin(trkhitx[i][ihitplane][j]);
            float Cx = X_correction_hists[ihitplane]->GetBinContent(x_bin);

            float Cyz = YZ_hist->GetBinContent(
                YZ_hist->FindBin(trkhitz[i][ihitplane][j],
                                                trkhity[i][ihitplane][j]));
	    
	    float Calibration_constant = CalibConstants[ihitplane]; // == SB debug, call calibration constant using plane number [ihitplane]
	    
            float corrected_dq_dx = trkdqdx[i][ihitplane][j]*Cx*norm_factors[ihitplane]*Cyz / Calibration_constant;
	    float dEdx_measured = corrected_dq_dx * Wion;
	    float total_field = tot_Ef(trkhitx[i][ihitplane][j],trkhity[i][ihitplane][j],trkhitz[i][ihitplane][j]);
	    /*
	    float dEdx_boxmodel = Dedx(corrected_dq_dx,
				       tot_Ef(trkhitx[i][ihitplane][j],
					      trkhity[i][ihitplane][j],
					      trkhitz[i][ihitplane][j]));
	    */

	    float KE_trk_hit = sp -> Eval(trkresrange[i][ihitplane][j]);
	    float dEdx_theory = dpdx(KE_trk_hit, trkpitch[i][ihitplane][j], Mmu);
	    
	    //cout << "SB debug corrected_dq_dx : " << corrected_dq_dx << ", dEdx_measured : " << dEdx_measured << ", dEdx_boxmodel : " << dEdx_boxmodel << ", dEdx_theory : " << dEdx_theory << ", total_field : " << total_field << endl;
	    //cout << "SB debug total_field : " << total_field << endl;
	    dedx_theory[ihitplane] -> Fill(dEdx_theory / total_field);
	    dedx_measured[ihitplane] -> Fill(dEdx_measured);
	    
	    if(bin_start > (dEdx_theory / total_field)) continue;
	    size_t bin = size_t(((dEdx_theory / total_field) - bin_start) / binsize);
	    //double bin_double = ((dEdx_theory / total_field) - bin_start) / binsize;

	    //cout << "SB debug bin : " << bin << ", bin_double : " << bin_double << endl;

	    if(bin < nbin){
	      dedx_measured_hist[ihitplane][bin] -> Fill(dEdx_measured);
	      dedx_theory_hist[ihitplane][bin] -> Fill(dEdx_theory);

	      dedx_measured_value[ihitplane][bin].push_back(dEdx_measured);
	      dedx_theory_value[ihitplane][bin].push_back(dEdx_theory);
	    }	    
          } // y containment.....
        } // loop over hits....
      }
    } // loop over crossing trks.......
  } // loop over jentries...........

  std::cout << "************************** Estimate dE/dx Medians and their error  *****************************" << std::endl;

  for(size_t ihitplane = 0; ihitplane < 3; ihitplane++) {
    for(size_t ibin = 0; ibin < nbin; ibin++){
      float this_theory_median = TMath::Median(dedx_theory_value[ihitplane][ibin].size(), &dedx_theory_value[ihitplane][ibin][0]);
      float this_measured_median = TMath::Median(dedx_measured_value[ihitplane][ibin].size(), &dedx_measured_value[ihitplane][ibin][0]);
      vector<float> current_theory_MAD;
      vector<float> current_measured_MAD;
      for(size_t ihit = 0; ihit < dedx_theory_value[ihitplane][ibin].size(); ihit++){
	current_theory_MAD.push_back(fabs(dedx_theory_value[ihitplane][ibin].at(ihit) - this_theory_median));
	current_measured_MAD.push_back(fabs(dedx_measured_value[ihitplane][ibin].at(ihit) - this_measured_median));
      }
      float this_theory_median_err = TMath::Median(current_theory_MAD.size(), &current_theory_MAD[0]);
      float this_measured_median_err = TMath::Median(current_measured_MAD.size(), &current_measured_MAD[0]);

      dedx_theory_median[ihitplane][ibin] = this_theory_median;
      dedx_measured_median[ihitplane][ibin] = this_measured_median;
      dedx_theory_median_err[ihitplane][ibin] = this_theory_median_err;
      dedx_measured_median_err[ihitplane][ibin] = this_measured_median_err;
      
      //cout << "SB debug, this_theory_median : " << this_theory_median << " +- " << this_theory_median_err << ", this_measured_median : " << this_measured_median << " +- " << this_measured_median_err << endl; 

    }
  }

  //vector<vector<float>> Recom_factor(3);
  //vector<vector<float>> Recom_factor_err(3);
  for(size_t ihitplane = 0; ihitplane < 3; ihitplane++) {
    for(size_t ibin = 0; ibin < nbin; ibin++){
      float this_Recom_fator = dedx_theory_median[ihitplane][ibin] / dedx_measured_median[ihitplane][ibin];
      float this_err = this_Recom_fator * sqrt(pow(dedx_theory_median_err[ihitplane][ibin]/dedx_theory_median[ihitplane][ibin], 2) + pow(dedx_measured_median_err[ihitplane][ibin]/dedx_measured_median[ihitplane][ibin], 2));
      Recom_factor[ihitplane][ibin] = this_Recom_fator;
      Recom_factor_err[ihitplane][ibin] = this_err;

      cout << "SB debug, (ihitplane, ibin) = (" << ihitplane << ", " << ibin << "), this_Recom_fator : " << this_Recom_fator << " +- " << this_err << endl;
    }
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  TH1D* this_template = new TH1D("", "", 1, 0., 4.);
  this_template -> GetYaxis() -> SetRangeUser(0.5, 1.5);
  this_template -> SetTitle("");
  this_template -> SetStats(0);
  this_template -> Draw();
  
  vector<float> x_axis;
  vector<float> x_axis_err;
  for(size_t ibin = 0; ibin < nbin; ibin++){
    float this_x = 3.025 + 0.05 * (ibin + 0.);
    x_axis.push_back(this_x);
    x_axis_err.push_back(0.025);
  }

  TGraphErrors *gr_plane0 = new TGraphErrors(nbin, &(x_axis[0]), &(Recom_factor[0][0]), &(x_axis_err[0]), &(Recom_factor_err[0][0]));
  TGraphErrors *gr_plane1 = new TGraphErrors(nbin, &(x_axis[0]), &(Recom_factor[1][0]), &(x_axis_err[0]), &(Recom_factor_err[1][0]));
  TGraphErrors *gr_plane2 = new TGraphErrors(nbin, &(x_axis[0]), &(Recom_factor[2][0]), &(x_axis_err[0]), &(Recom_factor_err[2][0]));

  gr_plane0 -> SetLineColor(kBlue);
  gr_plane1 -> SetLineColor(kRed);
  gr_plane2 -> SetLineColor(kGreen);

  gr_plane0 -> Draw("same");
  gr_plane1 -> Draw("same");
  gr_plane2 -> Draw("same");

  c -> SaveAs("../plots/Recom_fit_Birks.pdf");


  ////////////////////////////////////// Fitting Landau+Gaussian function to the histogram ////////////////////////////////
  //TSpline3 *sp = new TSpline3("Cubic Spline", &spline_Range[0], &spline_KE[0],13,"b2e2",0,0);

  outfile.cd();

  TVectorD NDF(3);
  TVectorD meanKE(3);
  TVectorD meanPitch(3);
  
  outfile.Write();
}

int main(int argc, char ** argv) {


  bool found_input = false,
       found_fcl = false,
       found_out = false;
  //string infile = argv[1];
  string infile;
  std::string fcl_file;

  double ref_dQdx[3] = {65.75, 63.5, 59.29};

  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
     found_fcl = true;
    }
    if (!strcasecmp(argv[iArg],"-i")) {
     infile = argv[++iArg];
     found_input = true;
    }
    if (!strcasecmp(argv[iArg],"-o")) {
      outfile_name = argv[++iArg];
      found_out = true;
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: dEdX_calibration " <<
                   "-i <input root file OR input list>" <<
                   "-c fclfile.fcl " << 
                   "-o outputfile.txt " << std::endl;
      return 1;
    }
  }

  if (!found_input) {
    cout << "Error: No input file was provided! Please provide with '-i'" << endl;
    return 0;
  }
  if (!found_fcl) {
    cout << "Error: No fcl file was provided! Please provide with '-c'" << endl;
    return 0;
  }
  if (!found_out) {
    cout << "Error: No output file was provided! Please provide with '-o'" << endl;
    return 0;
  }

  ////Setting up fcl parameters
  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;

  if (!fhicl_env) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing " <<
                 "or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
  auto const pset = fhicl::ParameterSet::make(fcl_file, lookupPolicy);
  /////

  //here
  std::string field_file = pset.get<std::string>("FieldMap");
  TFile * ef = new TFile("$DUNE_PARDATA_DIR/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root", "OPEN");
  //TFile * ef = OpenFile(field_file);

  std::vector<double> norm_factors = pset.get<std::vector<double>>("NormFactors");
  for (size_t i = 0; i<norm_factors.size(); ++i){
    if (norm_factors[i]>0){
      norm_factors[i] = ref_dQdx[i]/norm_factors[i];
    }
    else{
      norm_factors[i] = 1;
    }
    cout<<"Normalization factor["<<i<<"] = "<<norm_factors[i]<<endl;
  }
  std::vector<double> CalibConstants = pset.get<std::vector<double>>("CalibConstants");
  for (size_t i = 0; i < CalibConstants.size(); i++){
    cout<< "Calibration constant[" << i << "] = " << CalibConstants[i] << endl;
  }

  std::vector<std::pair<int, double>> KE_Range
      = pset.get<std::vector<std::pair<int, double>>>("KE_Range");
  //std::vector<int> spline_KE;
  //std::vector<double> spline_Range;
  for (size_t i = 0; i < KE_Range.size(); ++i) {
    spline_KE[i] = KE_Range[i].first;
    spline_Range[i] = KE_Range[i].second;
  }

  //usemap=true;

  TChain* shtree = new TChain("Event");
  TChain* shtree1 = new TChain("Event");
  TChain* shtree2 = new TChain("Event");

  if (infile.substr(infile.find_last_of(".") + 1) == "root"){
    shtree->Add(Form("%s/michelremoving2/Event", infile.c_str()));
    shtree1->Add(Form("%s/michelremoving2/Event", infile.c_str()));
    shtree2->Add(Form("%s/michelremoving2/Event", infile.c_str()));
  }
        
  else /*if(infile.substr(infile.find_last_of(".") + 1) == "txt")*/{
    std::ifstream in;
    in.open(infile.c_str());
    char line[1024];

    while(1){
     in.getline(line,1024);
      if (!in.good()) break;
      shtree->Add(Form("%s/michelremoving2/Event", line));
      shtree1->Add(Form("%s/michelremoving2/Event", line));
      shtree2->Add(Form("%s/michelremoving2/Event", line));
    }
    in.close();
    in.clear();
  }
  
  TFile output_file(outfile_name.c_str(), "RECREATE");

  //hitplane = 0;
  protoDUNE_Recom_fit *t=new protoDUNE_Recom_fit(shtree);
  t->GetEFMaps(ef);
  t->SetCaloMaps(pset);
  
  std::string method = pset.get<std::string>("Method","Lite");
  sce = pset.get<std::string>("SCE","on");
  if (sce=="on"){
    std::cout<<"SCE on"<<std::endl;
  }
  else{
    std::cout<<"SCE off"<<std::endl;
  }
  if (method=="Lite"){
    t->LoopLite(norm_factors,
                CalibConstants,
                output_file);
  }
  else{
    std::cout<<"Unknown method "<<method<<std::endl;
  }
  delete t;
  std::cout << "************************* Start of hitplane 1 ***************************" << std::endl;
  output_file.Close();
  
  /*
  //hitplane = 1;
  protoDUNE_dEdx_calib *t1=new protoDUNE_dEdx_calib(shtree1);
  t1->GetEFMaps(ef);
  t1->SetCaloMaps(pset);
  //for(calib_factor = 1.025e-3;  calib_factor<1.026e-3; calib_factor+=.001e-3) t1->Loop();
  for(double calib_factor = plane_1_low; calib_factor<plane_1_high; calib_factor+=plane_1_diff) t1->Loop(1, norm_factors[1], calib_factor);
  delete t1;
  std::cout << "************************* Start of hitplane 2 ***************************" << std::endl;

  //hitplane = 2;
  protoDUNE_dEdx_calib *t2=new protoDUNE_dEdx_calib(shtree2);
  t2->GetEFMaps(ef);
  t2->SetCaloMaps(pset);
  //for(calib_factor = 1.011e-3; calib_factor<1.012e-3; calib_factor+=.001e-3) t2->Loop();
  for(double calib_factor = plane_2_low; calib_factor<plane_2_high; calib_factor+=plane_2_diff) t2->Loop(2, norm_factors[2], calib_factor);
  delete t2; 
  */

} // main
