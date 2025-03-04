#ifndef PDSPTHINSLICEFITTER_hh
#define PDSPTHINSLICEFITTER_hh

#include <string>
#include <vector>
#include <map>

#include "TTree.h"
#include "TArrayD.h"
#include "TFile.h"
#include "THStack.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompChol.h"


#include "ThinSliceSample.h"
#include "ThinSliceDataSet.h"
#include "ThinSliceDriver.h"
#include "ThinSliceSystematic.h"
#include "ThinSliceEvent.h"

namespace protoana {

class PDSPThinSliceFitter {
 public:
  PDSPThinSliceFitter(std::string fcl_file, std::string output_file,
                      std::string mc_file = "", std::string data_file = "");
  void FillMCEvents();
  void BuildMCSamples();
  void SaveMCSamples();
  void GetNominalFluxes();
  void BuildDataHists();
  void InitializeMCSamples();
  void CompareDataMC(
      std::string extra_name, TDirectory * xsec_dir, TDirectory * plot_dir,
      bool post_fit = false);
  void ScaleMCToData();
  void ScaleDataToNorm();
  void RunFitAndSave();
  ~PDSPThinSliceFitter();

 private:
  void NormalFit();
  void Pulls();
  void Configure(std::string fcl_file);
  void DefineFitFunction();
  void MakeMinimizer();
  void ParameterScans();
  void DoThrows(const TH1D & pars, const TMatrixD * cov);
  void Do1DShifts(const TH1D & pars, bool prefit=false);
  void SetBestFit();
  void GetCurrentTruthHists(
    std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<TH1*>> & throw_inc_hists,
    std::map<int, std::vector<TH1*>> & throw_xsec_hists);
  void PlotThrows(
    std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<TH1*>> & truth_throw_hists,
    std::map<int, std::vector<TH1*>> & truth_inc_hists,
    std::map<int, std::vector<TH1*>> & truth_xsec_hists);
  void BuildFakeDataXSecs(bool use_scales = true);
  void BuildDataFromToy();
  double CalcChi2SystTerm(), CalcRegTerm();
  void MakeThrowsTree(TTree & tree, std::vector<double> & branches);
  //void MakeThrowsArrays(std::vector<TVectorD *> & arrays);

  std::vector<double> GetBestFitParsVec();

  ThinSliceDriver * fThinSliceDriver;
  std::map<int, std::vector<std::vector<ThinSliceSample>>> fSamples,
                                                           fFakeSamples;
  ThinSliceDataSet fDataSet;
  std::map<int, bool> fIsSignalSample;
  TFile fMCFile;
  //std::string fMCXSecFileName = "";
  TTree * fMCTree;
  TFile fDataFile;
  TTree * fDataTree;
  TFile fOutputFile;
  ROOT::Math::Functor fFitFunction;
  std::unique_ptr<ROOT::Math::Minimizer> fMinimizer;
  std::vector<double> fMinimizerInitVals;
  TH1D fPreFitParsNormal;

  //std::map<int, TGraphAsymmErrors> fSignalEfficiencies;
  //std::map<int, std::pair<TH1D, TH1D>> fSignalEffParts;

  //TGraphAsymmErrors fIncidentEfficiency;
  //TH1D fIncidentTotal, ;

  //std::map<int, TH1D> fSelectedDataHists;
  //std::map<int, TH1D> fRebinnedSelectedDataHists;
  //TH1D fIncidentDataHist;
  //TH1D fRebinnedIncidentDataHist;

  //THStack * fNominalIncidentMCStack;
  //THStack * fPostFitIncidentMCStack;
  //std::map<int, THStack *> fNominalSelectedMCStacks;
  //std::map<int, THStack *> fPostFitSelectedMCStacks;

  std::map<int, double> fNominalFluxes, fFakeFluxes;
  std::map<int, std::vector<std::vector<double>>> fFluxesBySample,
                                                  fFakeFluxesBySample;
  std::map<int, std::vector<int>> fFluxParsToSamples;
  double fDataFlux;
  double fMCDataScale = 1.;

  std::map<int, std::vector<double>> fSignalParameters;
  std::map<int, std::vector<std::string>> fSignalParameterNames;
  size_t fTotalSignalParameters = 0;

  std::map<int, double> fFluxParameters;
  std::map<int, std::string> fFluxParameterNames;
  size_t fTotalFluxParameters = 0;

  //std::map<int, std::string> fSystParameterNames;
  std::map<std::string, ThinSliceSystematic> fSystParameters;
  std::vector<std::string> fSystParameterNames;
  std::vector<double> fParLimits, fParLimitsUp;
  size_t fTotalSystParameters = 0;
  std::map<std::string, size_t> fCovarianceBins;
  bool fAddSystTerm, fAddRegTerm, fAddDiffInQuadrature;
  double fRegFactor = 0.;
  TMatrixD * fCovMatrix, * fCovMatrixDisplay;
  TDecompChol * fInputChol;
  //std::map<int, std::string> fDiffGraphs;
  std::string fDiffGraphFile, fDiffCovName;
  TH2D * fDiffCov;

  std::map<std::string, double> fToyValues;

  TRandom3 fRNG;
  std::map<int, std::vector<double>> fFakeDataScales;
  std::map<int, std::vector<double>> fBestFitSignalPars;
  std::map<std::string, ThinSliceSystematic> fBestFitSystPars;
  std::map<int, double> fBestFitFluxPars;
  std::map<int, TH1*> fNominalXSecs, fNominalIncs;
  std::map<int, TH1*> fBestFitXSecs, fBestFitIncs;
  std::map<int, TH1*> fFakeDataXSecs, fFakeDataIncs;

  std::map<int, TH1D*> fBestFitSelectionHists;
  std::map<int, std::vector<double>> fBestFitTruthVals;

  std::vector<ThinSliceEvent> fEvents, fFakeDataEvents;

  //Configurable members
  std::string fMCFileName;
  std::string fDataFileName;
  std::string fTreeName;
  std::vector<fhicl::ParameterSet> fSelectionSets;
  std::vector<fhicl::ParameterSet> fSampleSets;
  std::map<int, std::string> fFluxTypes;
  int fMaxCalls, fMaxIterations;
  size_t fNFitSteps = 0;
  unsigned int fNScanSteps;
  double fTolerance, fLowerLimit, fUpperLimit;
  std::vector<std::pair<int, int>> fPlotStyle;
  bool fPlotRebinned;
  bool fRandomStart;
  std::string fDriverName;
  std::string fAnalysis;
  fhicl::ParameterSet fAnalysisOptions;
  double fPitch;
  std::string fSliceMethod;
  bool fMultinomial;
  bool fDoFakeData, fDoThrows, fDoScans, fDo1DShifts, fDoSysts, fRunHesse;
  bool fFixVariables;
  std::map<std::string, double> fSystsToFix, fFixSystsPostFit;
  std::map<std::string, int> fSystParameterIndices;
  std::string fFakeDataRoutine;
  bool fDoFluctuateStats;
  bool fSplitMC;
  int fMaxEntries = -1, fMaxDataEntries = -1;
  int fSplitVal = 0;
  bool fFillIncidentInFunction = false;
  bool fFitUnderOverflow = false;
  bool fGetMeanXSec = false;
  bool fTieUnderOver = false;
  bool fUseFakeSamples = false;
  bool fFitFlux;
  size_t fNThrows, fMaxRethrows;
  std::string fFitType = "Normal";
  size_t fNPulls;
  bool fDoScaleDataToNorm;
  double fDataNorm;
  bool fDebugMCDataScale, fScaleToDataBeamProfile;
  bool fDebugChi2;

  std::vector<double> fIncidentRecoBins, fTrueIncidentBins, fBeamEnergyBins;
  std::vector<double> fDataBeamFluxes;
  std::vector<int> fIncidentSamples, fMeasurementSamples;
  bool fDrawXSecUnderflow;
  std::map<int, std::vector<double>> fSignalBins;
  //////////////////////////
  
  double BetheBloch(double energy, double mass) {
   double K,rho,Z,A, charge, me, I, gamma,  /*momentum ,*/wmax, pitch;
   long double beta;
    K = 0.307;
    rho = 1.4;
    charge = 1;
    Z = 18;
    A = 39.948;
    I = pow(10,-6)*10.5*18; //MeV
    me = 0.51; //MeV me*c^2
    pitch = 1;
    
    //momentum = sqrt( pow(energy,2) - pow(massicle,2));
    //beta = momentum/sqrt(pow(massicle,2) + pow(momentum,2));
    //gamma =  1/sqrt(1 - pow(beta,2));
    
    gamma = (energy + mass) / mass;
    beta = sqrt( 1 - 1/pow(gamma,2));

    wmax = 2*me*pow(beta,2)*pow(gamma,2)/(1+2*gamma*me/mass + pow(me,2)/pow(mass,2));
    
    
    double dEdX;
    //multiply by rho to have dEdX MeV/cm in LAr

    dEdX = pitch*(rho*K*Z*pow(charge,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma )/2 );

   return dEdX;
  };

  double densityEffect(double beta, double gamma) {
   double lar_C = 5.215, lar_x0 = 0.201, lar_x1 = 3, lar_a = 0.196, lar_k = 3;
   long double x = log10(beta * gamma);
   
   if( x >= lar_x1 ) return 2*log(10)*x - lar_C;

   else if ( lar_x0 <= x && x < lar_x1) return 2*log(10)*x - lar_C + lar_a * pow(( lar_x1 - x ) , lar_k );

   else return  0.; //if x < lar_x0
  };
};

}
#endif
