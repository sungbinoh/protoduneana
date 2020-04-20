////////////////////////////////////////////////////////////////////////
// Class:       PionAna
// Plugin Type: analyzer (art v2_07_03)
// File:        PionAna_module.cc
//
// Generated at Mon Sep  4 06:55:33 2017 by Leigh Whitehead using cetskelgen
// from cetlib version v3_00_01.
//
// This module is designed to show some usage examples of the analysis
// tools that I have been producing for protoDUNE. The aim has been to
// simplify the associations between the different objects to make
// some of the low-level art features more transparent
//
// The code is split into a few sections with each focusing on a different
// type of initial object
//
////////////////////////////////////////////////////////////////////////


#if defined(__MAKECINT__)
#pragma link C++ class vector<TVector3 >+; 
#endif

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
// #include "art/Framework/Services/Optional/TFileService.h"
#include "art_root_io/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "protoduneana/protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNECalibration.h"

 #include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"


// ROOT includes
#include "TTree.h"
#include "TVector3.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TLorentzVector.h" 


namespace protoana {
  class PionAna;
}


class protoana::PionAna : public art::EDAnalyzer {
public:

  explicit PionAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PionAna(PionAna const &) = delete;
  PionAna(PionAna &&) = delete;
  PionAna & operator = (PionAna const &) = delete;
  PionAna & operator = (PionAna &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  

  // Required functions.
  void analyze(art::Event const & e) override;




private:
  
  

  // fcl parameters
  std::string fParticleIDTag;
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fBeamModuleLabel;
  std::string fMCParticleTag;
  
  bool fVerbose;

  // struct BrokenTrack {
  //   const recob::Track * firstTrack;
  //   const recob::Track * secondTrack;
  //   double CosTheta;
  //   std::vector< float > Combined_ResidualRange;
  //   std::vector< float > Combined_dQdx;
  //   std::vector< float > Combined_dEdx;
  //   bool Valid;
  // };

//shared varibles

  int fEvent;        ///< number of the event being processed
  int fRun;          ///< number of the run being processed
  int fSubRun;       ///< number of the sub-run being processed
  int fNumDaughters; //
  int fNumBeamParticles;
  


  std::vector<float> fprimtrk_dqdx;  //removed outside vector as I only take one TPC track per event, need to emsure this is the case
  std::vector<float> fprimtrk_resrange;
  std::vector<float> fprimtrk_dedx;
  std::vector<float> fprimtrk_dedx_cal;
  // std::vector<float> fprimtrk_dedx_cal1;
  std::vector<TVector3> fprimtrk_hitpos;
  std::vector<float> fprimtrk_pitch;
  double fprimtrk_calorange; //track range from calo module. Does this match normal range?
  


  bool fIsBroken;
  double fCosTheta_brk; // if track ends in broken track region 

  std::vector< float > fCombined_ResidualRange;
  std::vector< float > fCombined_dQdx;
  std::vector< float > fCombined_dEdx;
  TVector3 fTrack_Start_brk;
  TVector3 fTrack_StartDir_brk;
  TVector3 fTrack_End_brk;
  TVector3 fTrack_EndDir_brk;
  float fTrack_pathlen_brk;
  float fTrack_calorange_brk;

  float fTrack_pathlen;
  TVector3 fTrack_Start;
  TVector3 fTrack_StartDir;
  TVector3 fTrack_End;
  TVector3 fTrack_EndDir;
  TVector3 fvtx;
  TVector3 finteractionvtx;



  int fBeamParticlePDG; //211 for track, 11 for shower
  char fBeamParticleType; //Track like or showerlike
  

  


  //Data only variobes

  std::vector<TVector3> fBeamEndPoints_BL;
  std::vector<TVector3> fBeamStartPoints_BL;
  std::vector<TVector3> fBeamDirectionStart_BL;


  int fNumParticles_BL; //BL -- beam line



  double fTOF_BL;
  std::vector<double> fMomentum_BL;
  std::vector< int > fPids;

  int fPID_Pdg[3];
  int fPID_Ndf[3];
  double fPID_MinChi2[3];
  double fPID_DeltaChi2[3];
  double fPID_Chi2Proton[3];
  double fPID_Chi2Kaon[3];
  double fPID_Chi2Pion[3];
  double fPID_Chi2Muon[3];
  double fPID_MissingE[3];
  double fPID_MissingEavg[3];
  double fPID_PIDA[3];     
//dsjd

  //MC only varibles
  int fBeamParticleTrackID_MCP;
  int fBeamParticlePDG_MCP;
  int fBeamParticleNumDaughters_MCP;
  int fBeamParticleMother_MCP;

  TVector3 fBeamParticleStart_MCP;
  TVector3 fBeamParticleStartMom_MCP;
  TVector3 fBeamParticleEnd_MCP;
  TVector3 fBeamParticleEndMom_MCP;

  std::string fBeamParticleEndProc_MCP;
  std::string fBeamParticleProc_MCP;


  std::vector<int> fDaughterParticleTrackID_MCP;
  std::vector<int> fDaughterParticlePDG_MCP;
  std::vector<int> fDaughterParticleNumDaughters_MCP;

  std::vector<TVector3> fDaughterParticleStart_MCP;
  std::vector<TVector3> fDaughterParticleStartMom_MCP;
  std::vector<TVector3> fDaughterParticleEnd_MCP;
  std::vector<TVector3> fDaughterParticleEndMom_MCP;
  std::vector<std::string> fDaughterParticleEndProc_MCP;
  std::vector<std::string> fDaughterParticleProc_MCP;


  bool fRecoTrackIsBeam_MCP;
  bool fFoundTrueMatch_MCP;




  int fBestRecoParticleTrackID_MCP;
  int fBestRecoParticlePDG_MCP;
  int fBestRecoParticleNumDaughters_MCP;
  int fBestRecoParticleMother_MCP;
  int fBestRecoParticleOrigin_MCP;
  TVector3 fBestRecoParticleStart_MCP;
  TVector3 fBestRecoParticleStartMom_MCP;
  TVector3 fBestRecoParticleEnd_MCP;
  TVector3 fBestRecoParticleEndMom_MCP;

  std::string fBestRecoParticleEndProc_MCP;
  std::string fBestRecoParticleProc_MCP;
  double fBestRecoParticleMatchFrac_MCP;

  std::vector<int> fRecoParticleTrackID_MCP;
  std::vector<int> fRecoParticlePDG_MCP;
  std::vector<int> fRecoParticleNumDaughters_MCP;

  std::vector<int> fRecoParticleOrigin_MCP;
  std::vector<int> fRecoParticleMother_MCP;

  std::vector<TVector3> fRecoParticleStart_MCP;
  std::vector<TVector3> fRecoParticleStartMom_MCP;
  std::vector<TVector3> fRecoParticleEnd_MCP;
  std::vector<TVector3> fRecoParticleEndMom_MCP;
  std::vector<std::string> fRecoParticleEndProc_MCP;
  std::vector<std::string> fRecoParticleProc_MCP;
  std::vector<double> fRecoParticleMatchFrac_MCP;



  TH1D* NumDaughters_hist;
  TH1D* TotalEvents_hist;
  TH1D* HasGoodBLInfo_hist;
  TH1D* HasGoodGeantPart_hist;
  TH1D* HasRecoTrack_hist;            
  TH1D* X_start;
  TH1D* NumBeamParticles;       

  TTree* EvtTree;     ///< tuple for each event
  TTree* CandidateTree; // tuple just for track like particles found

    //fhicl::ParameterSet beamlineUtil;
      
     //Beamline utils 
  protoana::ProtoDUNEBeamlineUtils beamlineUtil;
  protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNECalibration calibration;

  fhicl::ParameterSet fCalibrationPars;
  
  fhicl::ParameterSet fBrokenTrackParameters;

};


protoana::PionAna::PionAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fParticleIDTag(p.get<std::string>("ParticleIDTag")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fBeamModuleLabel(p.get<std::string>("BeamModuleLabel")),
  fMCParticleTag(p.get<std::string>("MCParticleTag")),
  fVerbose(p.get<bool>("Verbose")),
  beamlineUtil( p.get<fhicl::ParameterSet>("BeamlineUtils")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fCalibrationPars(p.get<fhicl::ParameterSet>("CalibrationPars")),
  
  fBrokenTrackParameters(p.get<fhicl::ParameterSet>("BrokenTrackParameters"))
{

}

void protoana::PionAna::beginJob()
{
  
  art::ServiceHandle<art::TFileService> tfs;
  NumDaughters_hist     = tfs->make<TH1D>("NumDaughters",";NumDaughters",  10, 0, 10);
  NumBeamParticles     = tfs->make<TH1D>("NumBeamParticles",";NumBeamParticles",10, 0, 10);

  TotalEvents_hist     = tfs->make<TH1D>("TotalEvents","TotalEvents;",3,0,3);
  HasGoodBLInfo_hist     = tfs->make<TH1D>("HasGoodBLInfo","HasGoodBLInfo;",3,0,3);
  HasGoodGeantPart_hist  = tfs->make<TH1D>("HasGoodGeantPart","HasGoodGeantPart;",3,0,3);
  HasRecoTrack_hist       = tfs->make<TH1D>("HasRecoTrack","HasRecoTrack;",3,0,3);
  
  // Define our n-tuples, which are limited forms of ROOT
  // TTrees. Start with the TTree itself.
  EvtTree    = tfs->make<TTree>("EvtTree",    "EvtTree");
  
  // Define the branches (columns) of our simulation n-tuple. To
  // write a variable, we give the address of the variable to
  // TTree::Branch.
  EvtTree->Branch("Event",       &fEvent,          "Event/I");
  EvtTree->Branch("SubRun",      &fSubRun,         "SubRun/I");
  EvtTree->Branch("Run",         &fRun,            "Run/I");
  EvtTree->Branch("NumBeamParticles",     &fNumBeamParticles,     "NumBeamParticles/I");

  EvtTree->Branch("NumDaughters",     &fNumDaughters,     "NumDaughters/I");
  EvtTree->Branch("BeamParticlePDG", &fBeamParticlePDG, "BeamParticlePDG/I");
  EvtTree->Branch("ParticleType",     &fBeamParticleType);
  EvtTree->Branch("NumParticles_BL", &fNumParticles_BL, "NumParticles_BL/I");
  EvtTree->Branch("vtx", &fvtx);
  EvtTree->Branch("TOF_BL", &fTOF_BL);
  EvtTree->Branch("Momentum_BL", &fMomentum_BL, 1,-1);
  EvtTree->Branch("Pids", &fPids);


  EvtTree->Branch("BeamParticleTrackID_MCP", &fBeamParticleTrackID_MCP, "BeamParticleTrackID_MCP/I");

  EvtTree->Branch("BeamParticlePDG_MCP", &fBeamParticlePDG_MCP, "BeamParticlePDG_MCP/I");
  EvtTree->Branch("BeamParticleStart_MCP", &fBeamParticleStart_MCP);
  EvtTree->Branch("BeamParticleStartMom_MCP", &fBeamParticleStartMom_MCP);
  EvtTree->Branch("BeamParticleEnd_MCP", &fBeamParticleEnd_MCP);
  EvtTree->Branch("BeamParticleEndMom_MCP", &fBeamParticleStartMom_MCP);
  EvtTree->Branch("BeamParticleNumDaughters_MCP", &fBeamParticleNumDaughters_MCP);
  EvtTree->Branch("BeamParticleEndProc_MCP", &fBeamParticleEndProc_MCP);
  EvtTree->Branch("BeamParticleProc_MCP", &fBeamParticleProc_MCP);

  
EvtTree->Branch("DaughterParticleTrackID_MCP", &fDaughterParticleTrackID_MCP,1,-1);
  EvtTree->Branch("DaughterParticlePDG_MCP", &fDaughterParticlePDG_MCP,1,-1);
  EvtTree->Branch("DaughterParticleStart_MCP", &fDaughterParticleStart_MCP,1,-1);
  EvtTree->Branch("DaughterParticleStartMom_MCP", &fDaughterParticleStartMom_MCP,1,-1);
  EvtTree->Branch("DaughterParticleEnd_MCP", &fDaughterParticleEnd_MCP,1,-1);
  EvtTree->Branch("DaughterParticleEndMom_MCP", &fDaughterParticleStartMom_MCP,1,-1);

EvtTree->Branch("DaughterParticleNumDaughters_MCP", &fDaughterParticleNumDaughters_MCP,1,-1);
EvtTree->Branch("DaughterParticleEndProc_MCP", &fDaughterParticleEndProc_MCP,1,-1);
EvtTree->Branch("DaughterParticleProc_MCP", &fDaughterParticleProc_MCP,1,-1);
  
  CandidateTree = tfs->make<TTree>("CandidateTree", "CandidateTree");

  
  CandidateTree->Branch("Event",       &fEvent,          "Event/I");
  CandidateTree->Branch("SubRun",      &fSubRun,         "SubRun/I");
  CandidateTree->Branch("Run",         &fRun,            "Run/I");

  CandidateTree->Branch("NumDaughters",     &fNumDaughters,     "NumDaughters/I");
  CandidateTree->Branch("BeamParticlePDG", &fBeamParticlePDG, "BeamParticlePDG/I");
  CandidateTree->Branch("ParticleType",     &fBeamParticleType);
  CandidateTree->Branch("Track_pathlen", &fTrack_pathlen, "Track_pathlen/F");
  
  CandidateTree->Branch("vtx", &fvtx);
  CandidateTree->Branch("interactionvtx", &finteractionvtx);
  CandidateTree->Branch("Track_Start", &fTrack_Start);
  CandidateTree->Branch("Track_End", &fTrack_End);
  CandidateTree->Branch("Track_StartDir", &fTrack_StartDir);
  CandidateTree->Branch("Track_EndDir", &fTrack_EndDir);

  CandidateTree->Branch("Track_pathlen_brk", &fTrack_pathlen_brk, "Track_pathlen_brk/F");
  CandidateTree->Branch("Track_calorange_brk", &fTrack_calorange_brk, "Track_calorange_brk/F");


  CandidateTree->Branch("Track_Start_brk", &fTrack_Start_brk);
  CandidateTree->Branch("Track_End_brk", &fTrack_End_brk);
  CandidateTree->Branch("Track_StartDir_brk", &fTrack_StartDir_brk);
  CandidateTree->Branch("Track_EndDir_brk", &fTrack_EndDir_brk);


  CandidateTree->Branch("IsBroken", &fIsBroken);
  
  CandidateTree->Branch("CosTheta_brk", &fCosTheta_brk);
  CandidateTree->Branch("Combined_ResidualRange", &fCombined_ResidualRange);
  CandidateTree->Branch("Combined_dQdx", &fCombined_dQdx);
  CandidateTree->Branch("Combined_dEdx", &fCombined_dEdx);



  CandidateTree->Branch("PID_Pdg",                &fPID_Pdg,               "PID_Pdg[3]/I");
  CandidateTree->Branch("PID_Ndf",                &fPID_Ndf,               "PID_Ndf[3]/I");
  CandidateTree->Branch("PID_MinChi2",            &fPID_MinChi2,           "PID_MinChi2[3]/D");
  CandidateTree->Branch("PID_DeltaChi2",          &fPID_DeltaChi2,         "PID_DeltaChi2[3]/D");
  CandidateTree->Branch("PID_Chi2Proton",         &fPID_Chi2Proton,        "PID_Chi2Proton[3]/D");
  CandidateTree->Branch("PID_Chi2Kaon",           &fPID_Chi2Kaon,          "PID_Chi2Kaon[3]/D");
  CandidateTree->Branch("PID_Chi2Pion",           &fPID_Chi2Pion,          "PID_Chi2Pion[3]/D");
  CandidateTree->Branch("PID_Chi2Muon",           &fPID_Chi2Muon,          "PID_Chi2Muon[3]/D");
  CandidateTree->Branch("PID_MissingE",           &fPID_MissingE,          "PID_MissingE[3]/D");
  CandidateTree->Branch("PID_MissingEavg",        &fPID_MissingEavg,       "PID_MissingEavg[3]/D");
  CandidateTree->Branch("PID_PIDA",               &fPID_PIDA,              "PID_PIDA[3]/D");  
  ///calo
  // CandidateTree->Branch("primtrk_dedx_cal1", &fprimtrk_dedx_cal1);
  CandidateTree->Branch("primtrk_dedx_cal", &fprimtrk_dedx_cal);
  CandidateTree->Branch("primtrk_dqdx", &fprimtrk_dqdx);
  CandidateTree->Branch("primtrk_dedx", &fprimtrk_dedx);
  CandidateTree->Branch("primtrk_calorange", &fprimtrk_calorange);
  CandidateTree->Branch("primtrk_resrange", &fprimtrk_resrange);
  CandidateTree->Branch("primtrk_pitch", &fprimtrk_pitch);
  CandidateTree->Branch("primtrk_hitpos", &fprimtrk_hitpos, 1,-1);
  
  CandidateTree->Branch("TOF_BL", &fTOF_BL);
  CandidateTree->Branch("Momentum_BL", &fMomentum_BL, 1,-1);
  CandidateTree->Branch("Pids", &fPids);
  CandidateTree->Branch("BeamEndPoints_BL", &fBeamEndPoints_BL, 1,-1);
  CandidateTree->Branch("BeamStartPoints_BL", &fBeamStartPoints_BL, 1,-1);
  CandidateTree->Branch("BeamDirectionStart_BL", &fBeamDirectionStart_BL, 1,-1);
  CandidateTree->Branch("NumParticles_BL", &fNumParticles_BL, "NumParticles_BL/I");
  
  CandidateTree->Branch("BeamParticleTrackID_MCP", &fBeamParticleTrackID_MCP, "BeamParticleTrackID_MCP/I");
  CandidateTree->Branch("BeamParticlePDG_MCP", &fBeamParticlePDG_MCP, "BeamParticlePDG_MCP/I");
  CandidateTree->Branch("BeamParticleStart_MCP", &fBeamParticleStart_MCP);
  CandidateTree->Branch("BeamParticleStartMom_MCP", &fBeamParticleStartMom_MCP);
  CandidateTree->Branch("BeamParticleEnd_MCP", &fBeamParticleEnd_MCP);
  CandidateTree->Branch("BeamParticleEndMom_MCP", &fBeamParticleStartMom_MCP);
  CandidateTree->Branch("BeamParticleNumDaughters_MCP", &fBeamParticleNumDaughters_MCP);
  CandidateTree->Branch("BeamParticleEndProc_MCP", &fBeamParticleEndProc_MCP);

  CandidateTree->Branch("BeamParticleProc_MCP", &fBeamParticleProc_MCP);
  CandidateTree->Branch("FoundTrueMatch_MCP", &fFoundTrueMatch_MCP);
  CandidateTree->Branch("RecoTrackIsBeam_MCP", &fRecoTrackIsBeam_MCP);

  CandidateTree->Branch("DaughterParticleTrackID_MCP", &fDaughterParticleTrackID_MCP,1,-1);
  CandidateTree->Branch("DaughterParticlePDG_MCP", &fDaughterParticlePDG_MCP,1,-1);
  CandidateTree->Branch("DaughterParticleStart_MCP", &fDaughterParticleStart_MCP,1,-1);
  CandidateTree->Branch("DaughterParticleStartMom_MCP", &fDaughterParticleStartMom_MCP,1,-1);
  CandidateTree->Branch("DaughterParticleEnd_MCP", &fDaughterParticleEnd_MCP,1,-1);
  CandidateTree->Branch("DaughterParticleEndMom_MCP", &fDaughterParticleStartMom_MCP,1,-1);
  CandidateTree->Branch("DaughterParticleNumDaughters_MCP", &fDaughterParticleNumDaughters_MCP,1,-1);
  CandidateTree->Branch("DaughterParticleEndProc_MCP", &fDaughterParticleEndProc_MCP,1,-1);
  CandidateTree->Branch("DaughterParticleProc_MCP", &fDaughterParticleProc_MCP,1,-1);

  CandidateTree->Branch("RecoParticleMother_MCP", &fRecoParticleMother_MCP,1,-1);
  CandidateTree->Branch("RecoParticleOrigin_MCP", &fRecoParticleOrigin_MCP,1,-1);
  CandidateTree->Branch("RecoParticleTrackID_MCP", &fRecoParticleTrackID_MCP,1,-1);
  CandidateTree->Branch("RecoParticlePDG_MCP", &fRecoParticlePDG_MCP,1,-1);
  CandidateTree->Branch("RecoParticleStart_MCP", &fRecoParticleStart_MCP,1,-1);
  CandidateTree->Branch("RecoParticleStartMom_MCP", &fRecoParticleStartMom_MCP,1,-1);
  CandidateTree->Branch("RecoParticleEnd_MCP", &fRecoParticleEnd_MCP,1,-1);
  CandidateTree->Branch("RecoParticleEndMom_MCP", &fRecoParticleEndMom_MCP,1,-1);
  CandidateTree->Branch("RecoParticleNumDaughters_MCP", &fRecoParticleNumDaughters_MCP,1,-1);
  CandidateTree->Branch("RecoParticleEndProc_MCP", &fRecoParticleEndProc_MCP,1,-1);
  CandidateTree->Branch("RecoParticleProc_MCP", &fRecoParticleProc_MCP,1,-1);
  CandidateTree->Branch("RecoParticleMatchFrac_MCP", &fRecoParticleMatchFrac_MCP,1,-1);


  CandidateTree->Branch("BestRecoParticleMother_MCP", &fBestRecoParticleMother_MCP);
  CandidateTree->Branch("BestRecoParticleOrigin_MCP", &fBestRecoParticleOrigin_MCP);
  CandidateTree->Branch("BestRecoParticleTrackID_MCP", &fBestRecoParticleTrackID_MCP);
  CandidateTree->Branch("BestRecoParticlePDG_MCP", &fBestRecoParticlePDG_MCP);
  CandidateTree->Branch("BestRecoParticleStart_MCP", &fBestRecoParticleStart_MCP);
  CandidateTree->Branch("BestRecoParticleStartMom_MCP", &fBestRecoParticleStartMom_MCP);
  CandidateTree->Branch("BestRecoParticleEnd_MCP", &fBestRecoParticleEnd_MCP);
  CandidateTree->Branch("BestRecoParticleEndMom_MCP", &fBestRecoParticleEndMom_MCP);
  CandidateTree->Branch("BestRecoParticleNumDaughters_MCP", &fBestRecoParticleNumDaughters_MCP);
  CandidateTree->Branch("BestRecoParticleEndProc_MCP", &fBestRecoParticleEndProc_MCP);
  CandidateTree->Branch("BestRecoParticleProc_MCP", &fBestRecoParticleProc_MCP);
  CandidateTree->Branch("BestRecoParticleMatchFrac_MCP", &fBestRecoParticleMatchFrac_MCP);
 // CandidateTree->Branch("BestAng", &fBestAng, "BestAng/F");


}




void protoana::PionAna::analyze(art::Event const & evt)
{

calibration = protoana::ProtoDUNECalibration( fCalibrationPars);

//clear at start as early returns mean end is not always reached

fBeamEndPoints_BL.clear();
fBeamDirectionStart_BL.clear();
fBeamStartPoints_BL.clear();
fMomentum_BL.clear();
fPids.clear();
fprimtrk_dedx_cal.clear();
fprimtrk_dqdx.clear();
fprimtrk_resrange.clear();
fprimtrk_dedx.clear();
fprimtrk_hitpos.clear();
fprimtrk_pitch.clear();



fDaughterParticleTrackID_MCP.clear();
fDaughterParticlePDG_MCP.clear();
fDaughterParticleNumDaughters_MCP.clear();
fDaughterParticleStart_MCP.clear();
fDaughterParticleStartMom_MCP.clear();
fDaughterParticleEnd_MCP.clear();
fDaughterParticleEndMom_MCP.clear();

fDaughterParticleEndProc_MCP.clear();
fDaughterParticleProc_MCP.clear();



fRecoParticleMother_MCP.clear();
fRecoParticleTrackID_MCP.clear();
fRecoParticlePDG_MCP.clear();
fRecoParticleNumDaughters_MCP.clear();
fRecoParticleStart_MCP.clear();
fRecoParticleStartMom_MCP.clear();
fRecoParticleEnd_MCP.clear();
fRecoParticleEndMom_MCP.clear();

fRecoParticleOrigin_MCP.clear();

fRecoParticleEndProc_MCP.clear();
fRecoParticleProc_MCP.clear();
fRecoParticleMatchFrac_MCP.clear();
//fTOF_BL.clear();


  fEvent  = evt.id().event(); 
  fRun    = evt.run();
  fSubRun = evt.subRun();
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  bool beamTriggerEvent = false;

  TotalEvents_hist->Fill(1);
  // If this event is MC then we can check what the true beam particle is
  if(!evt.isRealData()){
    // Get the truth utility to help us out
    protoana::ProtoDUNETruthUtils truthUtil;
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    // simulation has fGeneratorTag = "generator"
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
    // of these, so we pass the first element into the function to get the good particle
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);


    
    if(geantGoodParticle != 0x0){
      std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() << std::endl;
      HasGoodGeantPart_hist->Fill(1);
      fBeamParticleTrackID_MCP=geantGoodParticle->TrackId();
      fBeamParticlePDG_MCP=geantGoodParticle->PdgCode();
      fBeamParticleNumDaughters_MCP=geantGoodParticle->NumberDaughters();
      fBeamParticleStart_MCP=geantGoodParticle->Position(0).Vect();  //.Vect as i cant store TlorentzVector in std::vectpr
      fBeamParticleStartMom_MCP=geantGoodParticle->Momentum(0).Vect();
      fBeamParticleEnd_MCP=geantGoodParticle->EndPosition().Vect();
      fBeamParticleEndMom_MCP=geantGoodParticle->EndMomentum().Vect();
      fBeamParticleEndProc_MCP=geantGoodParticle->EndProcess();
      fBeamParticleProc_MCP=geantGoodParticle->Process();

      //largeant


      int nDaughters = geantGoodParticle->NumberDaughters();
      
      
      for(int iDaughter=0; iDaughter < nDaughters; iDaughter++){
        
        const int daughterID = geantGoodParticle->Daughter(iDaughter);
        
        
        
        const simb::MCParticle* Daughter = pi_serv->TrackIdToParticle_P(daughterID);
        //auto Daughter=getParticle(daughterID,evt);
        

        fDaughterParticleTrackID_MCP.push_back(Daughter->TrackId());
        fDaughterParticlePDG_MCP.push_back(Daughter->PdgCode());
        fDaughterParticleNumDaughters_MCP.push_back(Daughter->NumberDaughters());
        fDaughterParticleStart_MCP.push_back(Daughter->Position(0).Vect());
        fDaughterParticleStartMom_MCP.push_back(Daughter->Momentum(0).Vect());
        fDaughterParticleEnd_MCP.push_back(Daughter->EndPosition().Vect());
        fDaughterParticleEndMom_MCP.push_back(Daughter->EndMomentum().Vect());

        fDaughterParticleProc_MCP.push_back(Daughter->Process());
        fDaughterParticleEndProc_MCP.push_back(Daughter->EndProcess());
        
        }
      }//if geant particle exists
  }//is real data
  else{
    // For data we can see if this event comes from a beam trigger
    beamTriggerEvent = dataUtil.IsBeamTrigger(evt);
    if(beamTriggerEvent){
      std::cout << "This data event has a beam trigger" << std::endl;
    }



      //Access the Beam Event
  auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("beamevent");
  
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if( beamHandle.isValid()){
    art::fill_ptr_vector(beamVec, beamHandle);
  }

  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one
  /////////////////////////////////////////////////////////////
  
  
  //Check the quality of the event
  std::cout << "Timing Trigger: " << beamEvent.GetTimingTrigger() << std::endl; 
  std::cout << "Is Matched: "     << beamEvent.CheckIsMatched() << std::endl << std::endl;

  if( !beamlineUtil.IsGoodBeamlineTrigger(evt) ){
    std::cout << "Failed quality check" << std::endl;
    return;
  }
  HasGoodBLInfo_hist->Fill(1);

      std::vector<recob::Track> tracks= beamEvent.GetBeamTracks();
      
      fNumParticles_BL=tracks.size();
      
   


      for (size_t i = 0; i<tracks.size(); ++i){

        
        TVector3 pos_end = {tracks[i].End().X(), tracks[i].End().Y(), tracks[i].End().Z()}; 
        TVector3 pos_start= {tracks[i].Start().X(), tracks[i].Start().Y(), tracks[i].Start().Z()}; 
        TVector3 dir = {tracks[i].StartDirection().X(),tracks[i].StartDirection().Y(),tracks[i].StartDirection().Z()};
        fBeamEndPoints_BL.push_back(pos_end);
        fBeamStartPoints_BL.push_back(pos_start);
        fBeamDirectionStart_BL.push_back(dir);
      }



          //Get reconstructed beam momentum info
          auto & beammom = beamEvent.GetRecoBeamMomenta();
          for (size_t i = 0; i<beammom.size(); ++i){
            fMomentum_BL.push_back(beammom[i]);
      }

       //Access PID
  fPids = beamlineUtil.GetPID( beamEvent, 1. );


              if (beamEvent.GetTOFChan() != -1){//if TOFChan == -1, then there was not a successful match, if it's 0, 1, 2, or 3, then there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
            fTOF_BL = beamEvent.GetTOF();
        }
  }
  
  /*
  // Now we want to access the output from Pandora. This comes in the form of particle flow objects (recob::PFParticle).
  // The primary PFParticles are those we want to consider and these PFParticles then have a hierarchy of daughters that
  // describe the whole interaction of a given primary particle
  //
  //                     / daughter track
  //                    /
  //  primary track    /   
  //  ---------------- ---- daughter track
  //                   \
  //                   /\-
  //                   /\\-- daughter shower
  //
  // The above primary PFParticle will have links to three daughter particles, two track-like and one shower-like
  */

  // Get the PFParticle utility
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  fNumBeamParticles=beamParticles.size();
  NumBeamParticles->Fill(beamParticles.size());

  if(beamParticles.size() == 0){
    std::cerr << "We found no beam particles for this event... moving on" << std::endl;
    EvtTree->Fill();
    return;
  }
  
  
  // We can now look at these particles
  for(const recob::PFParticle* particle : beamParticles){

    // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
    // of this particle might be more helpful. These return null pointers if not track-like / shower-like
    fvtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fBeamParticlePDG=particle->PdgCode();
    const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
    
    if(thisShower != 0x0){
      std::cout << "Beam particle is shower-like, returning" << std::endl;
      fBeamParticleType='s';
      EvtTree->Fill();
      return;
    } 

    if(thisTrack == 0x0){
      std::cout << "No track, return" << std::endl;
      EvtTree->Fill();
      return;
    }

    else{
      std::cout << "Beam particle is track-like" << std::endl;
      fBeamParticleType='t';
    }





    // beam information
    // std::vector<TVector3> beamEndPos;
    // std::vector<TVector3> beamDir;

    HasRecoTrack_hist->Fill(1);
    
      
      // Find the particle vertex. We need the tracker tag here because we need to do a bit of
    // additional work if the PFParticle is track-like to find the vertex. 
    fvtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    
    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    finteractionvtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);

    

    // Get the track utility
    protoana::ProtoDUNETrackUtils trackUtil;

    TVector3 track_vector= finteractionvtx-fvtx;
    std::vector<anab::ParticleID> pids = trackUtil.GetRecoTrackPID(*thisTrack, evt, fTrackerTag, fParticleIDTag);

    
    for(size_t k = 0; k < pids.size() && k<3; k++){
      int plane = pids[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fPID_Pdg[plane]            = pids[plane].Pdg();
      fPID_Ndf[plane]            = pids[plane].Ndf();
      fPID_MinChi2[plane]        = pids[plane].MinChi2();
      fPID_DeltaChi2[plane]      = pids[plane].DeltaChi2();
      fPID_Chi2Proton[plane]     = pids[plane].Chi2Proton();
      fPID_Chi2Kaon[plane]       = pids[plane].Chi2Kaon();
      fPID_Chi2Pion[plane]       = pids[plane].Chi2Pion();
      fPID_Chi2Muon[plane]       = pids[plane].Chi2Muon();
      fPID_MissingE[plane]       = pids[plane].MissingE();
      fPID_MissingEavg[plane]    = pids[plane].MissingEavg();
      fPID_PIDA[plane]           = pids[plane].PIDA();
    }

    if (!evt.isRealData()){
          protoana::ProtoDUNETruthUtils truthUtil;
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    // simulation has fGeneratorTag = "generator"
          auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
          
    // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
    // of these, so we pass the first element into the function to get the good particle
          fRecoTrackIsBeam_MCP=false;
          fFoundTrueMatch_MCP=false;
          
        //   std::cout<<"origin1: "<<mcTruths[0].Origin()<<std::endl;
              
  
        



      
        // std::cout<<"origin2: "<<mcTruths[1].Origin()<<std::endl;
          
          
          const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
          auto trueTrack_vec=truthUtil.GetMCParticleListFromRecoTrack(*thisTrack,evt,fTrackerTag);
           //gets vector of pair of MCpartiles and their frac contribution 
          
          double BestMaxFrac=0;
          const simb::MCParticle* trueTrack = 0x0;
         
          for(size_t i = 0; i<trueTrack_vec.size(); ++i)
          {
          //or(auto pair : trueTrack_vec){
            
            fRecoParticlePDG_MCP.push_back(trueTrack_vec[i].first->PdgCode());
            // std::cout<<trueTrack_vec[i].first::std<<std::endl;
            
            fRecoParticleTrackID_MCP.push_back(trueTrack_vec[i].first->TrackId());
            
            
            fRecoParticleNumDaughters_MCP.push_back(trueTrack_vec[i].first->NumberDaughters());
            fRecoParticleStart_MCP.push_back(trueTrack_vec[i].first->Position(0).Vect());
            fRecoParticleStartMom_MCP.push_back(trueTrack_vec[i].first->Momentum(0).Vect());
            fRecoParticleEnd_MCP.push_back(trueTrack_vec[i].first->EndPosition().Vect());
            fRecoParticleEndMom_MCP.push_back(trueTrack_vec[i].first->EndMomentum().Vect());
            fRecoParticleProc_MCP.push_back(trueTrack_vec[i].first->Process());
            fRecoParticleEndProc_MCP.push_back(trueTrack_vec[i].first->EndProcess());
            fRecoParticleMother_MCP.push_back(trueTrack_vec[i].first->Mother());
            fRecoParticleMatchFrac_MCP.push_back(trueTrack_vec[i].second);

            
            auto mctruth=pi_serv->ParticleToMCTruth_P(trueTrack_vec[i].first);
            fRecoParticleOrigin_MCP.push_back(mctruth->Origin()); //take first mctut


            
            if (trueTrack_vec[i].second > BestMaxFrac){

              fBestRecoParticleMatchFrac_MCP=trueTrack_vec[i].second;
              BestMaxFrac=fBestRecoParticleMatchFrac_MCP;
              trueTrack=trueTrack_vec[i].first;
              fBestRecoParticleTrackID_MCP=trueTrack->TrackId();  
              fBestRecoParticlePDG_MCP=trueTrack->PdgCode();
              fBestRecoParticleNumDaughters_MCP=trueTrack->NumberDaughters();
              fBestRecoParticleStart_MCP=trueTrack->Position(0).Vect();
              fBestRecoParticleStartMom_MCP=trueTrack->Momentum(0).Vect();
              
              fBestRecoParticleEnd_MCP=trueTrack->EndPosition().Vect();
              
              fBestRecoParticleEndMom_MCP=trueTrack->EndMomentum().Vect();
              
              fBestRecoParticleProc_MCP=trueTrack->Process();
              
              fBestRecoParticleEndProc_MCP=trueTrack->EndProcess();
              
              fBestRecoParticleMother_MCP=trueTrack->Mother();
              
              auto mctruth=pi_serv->ParticleToMCTruth_P(trueTrack);
              
              fBestRecoParticleOrigin_MCP=mctruth->Origin();
              
              fFoundTrueMatch_MCP=true;
            }

         }
          
          if((geantGoodParticle != 0x0) & (trueTrack!= 0x0)){
          
          if(trueTrack->TrackId()==geantGoodParticle->TrackId()) fRecoTrackIsBeam_MCP=true;
          
        }// if both particles exits

    } //if not real data 

    
    
    fTrack_pathlen=track_vector.Mag();
    fTrack_Start = {thisTrack->Start().X(),thisTrack->Start().Y(),thisTrack->Start().Z()};
    fTrack_StartDir = {thisTrack->StartDirection().X(),thisTrack->StartDirection().Y(),thisTrack->StartDirection().Z()};
    fTrack_End = {thisTrack->End().X(),thisTrack->End().Y(),thisTrack->End().Z()};
    fTrack_EndDir = {thisTrack->EndDirection().X(),thisTrack->EndDirection().Y(),thisTrack->EndDirection().Z()};

    
 
    /// Try to determine if it's a broken track
    auto BrokenTrack = trackUtil.IsBrokenTrack(*thisTrack, evt, fTrackerTag,fCalorimetryTag,fBrokenTrackParameters,fCalibrationPars);
  
  fIsBroken=BrokenTrack.Valid;
  std::cout<<"Broken? "<<fIsBroken<<std::endl;
  fCosTheta_brk=BrokenTrack.CosTheta; // if track ends in broken track region 
  fCombined_ResidualRange=BrokenTrack.Combined_ResidualRange; 
  fCombined_dEdx=BrokenTrack.Combined_dEdx;
  fCombined_dQdx=BrokenTrack.Combined_dQdx;
 
  if(BrokenTrack.Valid){
    fTrack_pathlen_brk=BrokenTrack.secondTrack->Length();
    
    fTrack_Start_brk = {BrokenTrack.secondTrack->Start().X(),BrokenTrack.secondTrack->Start().Y(),BrokenTrack.secondTrack->Start().Z()};
    
    fTrack_StartDir_brk = {BrokenTrack.secondTrack->StartDirection().X(),BrokenTrack.secondTrack->StartDirection().Y(),BrokenTrack.secondTrack->StartDirection().Z()};
    
    fTrack_End_brk = {BrokenTrack.secondTrack->End().X(),thisTrack->End().Y(),BrokenTrack.secondTrack->End().Z()};
    fTrack_EndDir_brk = {BrokenTrack.secondTrack->EndDirection().X(),BrokenTrack.secondTrack->EndDirection().Y(),BrokenTrack.secondTrack->EndDirection().Z()};
    std::vector<anab::Calorimetry> calovector_brk = trackUtil.GetRecoTrackCalorimetry(*BrokenTrack.secondTrack, evt, fTrackerTag, fCalorimetryTag);
    
    for (auto & calo : calovector_brk){
       if (calo.PlaneID().Plane == 2){ //only collection plane
          fTrack_calorange_brk=calo.Range();
        }
      }

    
}


  std::vector<anab::ParticleID> daughterpids = trackUtil.GetRecoTrackPID(*thisTrack, evt, fTrackerTag, fParticleIDTag);
   


                    //HY::Get the Calorimetry(s) from thisTrack
    std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
 
    
    
    fprimtrk_dedx_cal = calibration.GetCalibratedCalorimetry(  *thisTrack, evt, fTrackerTag, fCalorimetryTag );
    // fprimtrk_dedx_cal=trackUtil.CalibrateCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag,fCalorimetryParameters);
    

    

    for (auto & calo : calovector){
       if (calo.PlaneID().Plane == 2){ //only collection plane
          fprimtrk_calorange=calo.Range();
          for (size_t ihit = 0; ihit < calo.dQdx().size(); ++ihit){ //loop over hits
              fprimtrk_dqdx.push_back(calo.dQdx()[ihit]);
              fprimtrk_resrange.push_back(calo.ResidualRange()[ihit]);
              fprimtrk_dedx.push_back(calo.dEdx()[ihit]);
              fprimtrk_pitch.push_back(calo.TrkPitchVec()[ihit]);
        TVector3 primtrk_hitpos_tmp;
        const auto &primtrk_pos=(calo.XYZ())[ihit];
              primtrk_hitpos_tmp.SetX(primtrk_pos.X()); //convert from stupid ROOT math 3d vector
              primtrk_hitpos_tmp.SetY(primtrk_pos.Y());
              primtrk_hitpos_tmp.SetZ(primtrk_pos.Z());
              fprimtrk_hitpos.push_back(primtrk_hitpos_tmp);

                            } //loop over hits
                       } //only collection plane
                }//loop over calo vector
           
    // Let's get the daughter PFParticles... we can do this simply without the utility
    for(const int daughterID : particle->Daughters()){
      // Daughter ID is the element of the original recoParticle vector
      const recob::PFParticle *daughterParticle = &(recoParticles->at(daughterID));
      std::cout << "Daughter " << daughterID << " has " << daughterParticle->NumDaughters() << " daughters" << std::endl;
    }
 
    // For actually studying the objects, it is easier to have the daughters in their track and shower forms.
    // We can use the utility to get a vector of track-like and a vector of shower-like daughters
    const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
    const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  
    std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;
   
    fNumDaughters=trackDaughters.size()+showerDaughters.size(); //this will get written multiple times in event of multiple beam particles
  }




NumDaughters_hist->Fill(fNumDaughters);
EvtTree->Fill();
CandidateTree->Fill();


}

void protoana::PionAna::endJob()
{

}

DEFINE_ART_MODULE(protoana::PionAna)

