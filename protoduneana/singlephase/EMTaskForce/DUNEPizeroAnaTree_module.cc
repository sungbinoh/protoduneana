////////////////////////////////////////////////////////////////////////
// Class:       DUNEPizeroAnaTree
// File:        DUNEPizeroAnaTree_module.cc
//
// Extract DUNE useful information, do a first pre-selection and
// save output to a flat tree
//
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
// modified by Aaron Higuera ahiguera@central.uh.edu
// and Milo Vermeulen milov@nikhef.nl
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcore/Geometry/Geometry.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dune/Calib/XYZCalib.h"
#include "dune/CalibServices/XYZCalibService.h"
#include "dune/CalibServices/XYZCalibServiceProtoDUNE.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "PiZeroProcess.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"

// C++ Includes
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_map>

// Maximum number of beam particle daughters to save
const int NMAXDAUGHTERS = 25;
// Maximum number of primary MC pi0s to save
const int NMAXPIZEROS = 25;
// Maximum number of hits to save
const int NMAXHITS = 5000;

namespace protoana {
  class DUNEPizeroAnaTree;
}


class protoana::DUNEPizeroAnaTree : public art::EDAnalyzer {
public:

  explicit DUNEPizeroAnaTree(fhicl::ParameterSet const & p);

  DUNEPizeroAnaTree(DUNEPizeroAnaTree const &) = delete;
  DUNEPizeroAnaTree(DUNEPizeroAnaTree &&) = delete;
  DUNEPizeroAnaTree & operator = (DUNEPizeroAnaTree const &) = delete;
  DUNEPizeroAnaTree & operator = (DUNEPizeroAnaTree &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & evt) override;

private:

  // Helper utility functions
  // protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNEShowerUtils showerUtil;
  // protoana::ProtoDUNEBeamlineUtils beamlineUtil;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;


  // Track momentum algorithm calculates momentum based on track range
  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // Initialise tree variables
  void Initialise();
  // void ResetPi0Vars();
  void FillPrimaryPFParticle(art::Event const & evt, const recob::PFParticle* particle);
  void FillAPA3Object(art::Event const & evt, art::Handle<std::vector<recob::PFParticle>> pfpHandle);
  void setPiZeroInfo(art::Event const & evt, const std::vector<const simb::MCParticle*>& pi0s);

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fSimulationTag;
  std::string fCalorimetryTag;
  calo::CalorimetryAlg fCalorimetryAlg;
  std::string fParticleIDTag;
  std::string fTrackTag;
  std::string fShowerTag;
  // std::string fShowerCaloTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fHitTag;
  int fVerbose;

  TTree* fNuBeam;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
  int fNactivefembs[6];

  // Beam track
  int fbeamtrigger;
  int fbeamCheckIsMatched;
  double ftof;
  int fcerenkovStatus[2];
  double fcerenkovTime[2];
  double fcerenkovPressure[2];
  double fbeamtrackMomentum;
  double fbeamtrackEnergy;
  double fbeamtrackPos[3];
  double fbeamtrackEndPos[3];
  double fbeamtrackDir[3];
  double fbeamtrackTime;
  int fbeam_ntrjPoints;
  int fbeamtrackPdg;
  int fbeamtrackID;
  double fbeamtrackPos_at[1000][3];
  double fbeamtrackMom_at[1000][4];
  std::vector<std::string> fbeamtrackEndProcess;
  int fbeamtrackNDaughters;

  // Reconstructed primary particle
  double fprimaryVertex[3];
  int fprimaryIstrack;
  int fprimaryIsshower;
  double fprimaryBDTScore;
  double fprimaryCNNScore;
  int fprimaryNHits;
  double fprimaryTheta;
  double fprimaryPhi;
  double fprimaryLength;
  double fprimaryMomentum;
  double fprimaryEndMomentum;
  double fprimaryEndPosition[3];
  double fprimaryStartPosition[3];
  double fprimaryEndDirection[3];
  double fprimaryStartDirection[3];
  double fprimaryOpeningAngle;
  int fprimaryShowerBestPlane;
  double fprimaryShowerEnergy;
  double fprimaryShowerCharge;
  double fprimaryShowerMIPEnergy;
  double fprimaryMomentumByRangeProton;
  int fprimaryIsBeamparticle;
  int    fprimaryTruth_trkID;
  int    fprimaryTruth_pdg;
  double fprimaryTruth_E;
  double fprimaryTruth_vtx[3];
  double fprimaryKineticEnergy[3];
  double fprimaryRange[3];
  int    fprimarynCal;
  double fprimarydEdx[NMAXHITS];
  double fprimarydQdx[NMAXHITS];
  double fprimary_cal_pos[NMAXHITS][3];
  double fprimary_cal_pitch[NMAXHITS];
  double fprimaryResidualRange[NMAXHITS];
  int    fprimaryShower_nHits; //collection only
  int    fprimaryShower_hit_w[NMAXHITS];
  double fprimaryShower_hit_q[NMAXHITS];
  double fprimaryShower_hit_t[NMAXHITS];
  double fprimaryShower_hit_pos[NMAXHITS][3];
  int fprimaryID;
  double fprimaryT0;
  int fprimaryNDaughters;

  // Primary particle daughters
  int fprimaryDaughterID[NMAXDAUGHTERS];
  int fprimaryDaughterIstrack[NMAXDAUGHTERS];
  int fprimaryDaughterIsshower[NMAXDAUGHTERS];
  double fprimaryDaughterBDTScore[NMAXDAUGHTERS];
  double fprimaryDaughterCNNScore[NMAXDAUGHTERS];
  int fprimaryDaughterNCollHits[NMAXDAUGHTERS];
  int fprimaryDaughterNHits[NMAXDAUGHTERS];
  double fprimaryDaughterEnergy[NMAXDAUGHTERS];
  double fprimaryDaughterEnergyFromHits[NMAXDAUGHTERS];
  double fprimaryDaughterMomentum[NMAXDAUGHTERS];
  double fprimaryDaughterEndMomentum[NMAXDAUGHTERS];
  double fprimaryDaughterLength[NMAXDAUGHTERS];
  double fprimaryDaughterStartPosition[NMAXDAUGHTERS][3];
  double fprimaryDaughterEndPosition[NMAXDAUGHTERS][3];
  double fprimaryDaughterStartDirection[NMAXDAUGHTERS][3];
  double fprimaryDaughterEndDirection[NMAXDAUGHTERS][3];
  int fprimaryDaughterParentPdg[NMAXDAUGHTERS];
  int fprimaryDaughterGrandparentPdg[NMAXDAUGHTERS];
  int fprimaryDaughterParentID[NMAXDAUGHTERS];
  int fprimaryDaughterGrandparentID[NMAXDAUGHTERS];
  double fprimaryDaughterParentE[NMAXDAUGHTERS];
  double fprimaryDaughterGrandparentE[NMAXDAUGHTERS];
  double fprimaryDaughterParentStart[NMAXDAUGHTERS][3];
  double fprimaryDaughterGrandparentStart[NMAXDAUGHTERS][3];
  double fprimaryDaughterParentEnd[NMAXDAUGHTERS][3];
  double fprimaryDaughterGrandparentEnd[NMAXDAUGHTERS][3];
  int fprimaryDaughterParentEndProcess[NMAXDAUGHTERS];
  double fprimaryDaughterCompleteness[NMAXDAUGHTERS];
  double fprimaryDaughterPurity[NMAXDAUGHTERS];
  // double fprimaryDaughterHitdQdx[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterCaloHitdQdx[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterCaloHitPitch[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterOwnHitPitch[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterHitCharge[NMAXDAUGHTERS][NMAXHITS];
  // double fprimaryDaughterHitdEdxArea[NMAXDAUGHTERS][NMAXHITS];
  int fprimaryDaughterHitWire[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterHitTime[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterHitPosition[NMAXDAUGHTERS][NMAXHITS][3];
  // double fprimaryDaughterCaloHitPosition[NMAXDAUGHTERS][NMAXHITS][3];

  // MC pi0 variables
  int fNMCPi0s;
  int fMCPi0ID[NMAXPIZEROS];
  double fMCPi0Energy[NMAXPIZEROS];
  double fMCPi0StartPosition[NMAXPIZEROS][3];
  double fMCPi0EndPosition[NMAXPIZEROS][3];
  int fMCPhoton1ID[NMAXPIZEROS];
  int fMCPhoton2ID[NMAXPIZEROS];
  double fMCPhoton1Energy[NMAXPIZEROS];
  double fMCPhoton2Energy[NMAXPIZEROS];
  double fMCPhoton1StartPosition[NMAXPIZEROS][3];
  double fMCPhoton2StartPosition[NMAXPIZEROS][3];
  double fMCPhoton1EndPosition[NMAXPIZEROS][3];
  double fMCPhoton2EndPosition[NMAXPIZEROS][3];
  // // Reco hits related to photons
  // int fMCPhoton1NumHits[NMAXPIZEROS];
  // int fMCPhoton2NumHits[NMAXPIZEROS];
  // int fMCPhoton1_hit_w[NMAXPIZEROS][NMAXHITS];
  // int fMCPhoton2_hit_w[NMAXPIZEROS][NMAXHITS];
  // double fMCPhoton1_hit_t[NMAXPIZEROS][NMAXHITS];
  // double fMCPhoton2_hit_t[NMAXPIZEROS][NMAXHITS];
  // double fMCPhoton1_hit_q[NMAXPIZEROS][NMAXHITS];
  // double fMCPhoton2_hit_q[NMAXPIZEROS][NMAXHITS];
  // double fMCPhoton1_hit_pos[NMAXPIZEROS][NMAXHITS][3];
  // double fMCPhoton2_hit_pos[NMAXPIZEROS][NMAXHITS][3];

  // Other
  geo::GeometryCore const* fGeometryService;
  spacecharge::SpaceChargeService::provider_type const* fSCEService;
};


protoana::DUNEPizeroAnaTree::DUNEPizeroAnaTree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  // dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  // beamlineUtil(p.get<fhicl::ParameterSet>("BeamLineUtils")),
  // fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fSimulationTag(p.get<std::string>("SimulationTag")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fParticleIDTag(p.get<std::string>("ParticleIDTag")),
  fTrackTag(p.get<std::string>("TrackTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  // fShowerCaloTag(p.get<std::string>("ShowerCalorimetryTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fVerbose(p.get<int>("Verbose")) {
  // Get a pointer to the geometry service provider.
  fGeometryService = lar::providerFrom<geo::Geometry>();
  fSCEService = lar::providerFrom<spacecharge::SpaceChargeService>();
}

void protoana::DUNEPizeroAnaTree::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fNuBeam = tfs->make<TTree>("PandoraBeam", "Beam events reconstructed with Pandora");
  fNuBeam->Branch("run",                           &fRun,                          "run/I");
  fNuBeam->Branch("subrun",                        &fSubRun,                       "subrun/I");
  fNuBeam->Branch("event",                         &fevent,                        "event/I");
  fNuBeam->Branch("timestamp",                     &fTimeStamp,                    "timestamp/D");
  fNuBeam->Branch("Nactivefembs",                  &fNactivefembs,                 "Nactivefembs[5]/I");

  fNuBeam->Branch("beamtrigger",                   &fbeamtrigger,                  "beamtrigger/I");
  fNuBeam->Branch("beamCheckIsMatched",            &fbeamCheckIsMatched,           "beamCheckIsMatched/I");
  fNuBeam->Branch("tof",                           &ftof,                          "tof/D");
  fNuBeam->Branch("cerenkovStatus",                &fcerenkovStatus,               "cerenkovStatus[2]/I");
  fNuBeam->Branch("cerenkovTime",                  &fcerenkovTime,                 "cerenkovTime[2]/D");
  fNuBeam->Branch("cerenkovPressure",              &fcerenkovPressure,             "cerenkovPressure[2]/D");
  fNuBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum/D");
  fNuBeam->Branch("beamtrackEnergy",               &fbeamtrackEnergy,              "beamtrackEnergy/D");
  fNuBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[3]/D");
  fNuBeam->Branch("beamtrackEndPos",               &fbeamtrackEndPos,              "beamtrackEndPos[3]/D");
  fNuBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[3]/D");
  fNuBeam->Branch("beamtrackTime",                 &fbeamtrackTime,                "beamtrackTime/D");
  fNuBeam->Branch("beamtrackPdg",                  &fbeamtrackPdg,                 "beamtrackPdg/I");
  fNuBeam->Branch("beamtrackID",                   &fbeamtrackID,                  "beamtrackID/I");
  fNuBeam->Branch("beam_ntrjPoints",               &fbeam_ntrjPoints,              "beam_ntrjPoints/I");
  fNuBeam->Branch("beamtrackNDaughters",           &fbeamtrackNDaughters,          "beamtrackNDaughters/I");
  fNuBeam->Branch("beamtrackPos_at",               &fbeamtrackPos_at,              "beamtrackPos_at[beam_ntrjPoints][3]/D");
  fNuBeam->Branch("beamtrackMom_at",               &fbeamtrackMom_at,              "beamtrackMom_at[beam_ntrjPoints][4]/D");
  fNuBeam->Branch("beamtrackEndProcess",           &fbeamtrackEndProcess);

  fNuBeam->Branch("primaryVertex",                 &fprimaryVertex,                "primaryVertex[3]/D");
  fNuBeam->Branch("primaryIsBeamparticle",         &fprimaryIsBeamparticle,        "primaryIsBeamparticle/I");
  fNuBeam->Branch("primaryIstrack",                &fprimaryIstrack,               "primaryIstrack/I");
  fNuBeam->Branch("primaryIsshower",               &fprimaryIsshower,              "primaryIsshower/I");
  fNuBeam->Branch("primaryBDTScore",               &fprimaryBDTScore,              "primaryBDTScore/D");
  fNuBeam->Branch("primaryCNNScore",               &fprimaryCNNScore,              "primaryCNNScore/D");
  fNuBeam->Branch("primaryNHits",                  &fprimaryNHits,                 "primaryNHits/I");
  fNuBeam->Branch("primaryTheta",                  &fprimaryTheta,                 "primaryTheta/D");
  fNuBeam->Branch("primaryPhi",                    &fprimaryPhi,                   "primaryPhi/D");
  fNuBeam->Branch("primaryLength",                 &fprimaryLength,                "primaryLength/D");
  fNuBeam->Branch("primaryMomentum",               &fprimaryMomentum,              "primaryMomentum/D");
  fNuBeam->Branch("primaryEndMomentum",            &fprimaryEndMomentum,           "primaryEndMomentum/D");
  fNuBeam->Branch("primaryEndPosition",            &fprimaryEndPosition,           "primaryEndPosition[3]/D");
  fNuBeam->Branch("primaryStartPosition",          &fprimaryStartPosition,         "primaryStartPosition[3]/D");
  fNuBeam->Branch("primaryEndDirection",           &fprimaryEndDirection,          "primaryEndDirection[3]/D");
  fNuBeam->Branch("primaryStartDirection",         &fprimaryStartDirection,        "primaryStartDirection[3]/D");
  fNuBeam->Branch("primaryOpeningAngle",           &fprimaryOpeningAngle,          "primaryOpeningAngle/D");
  fNuBeam->Branch("primaryID",                     &fprimaryID,                    "primaryID/I");
  fNuBeam->Branch("primaryTruth_E",                &fprimaryTruth_E,               "primaryTruth_E/D");
  fNuBeam->Branch("primaryTruth_vtx",              &fprimaryTruth_vtx,             "primaryTruth_vtx[3]/D");
  fNuBeam->Branch("primaryTruth_pdg",              &fprimaryTruth_pdg,             "primaryTruth_pdg/I");
  fNuBeam->Branch("primaryTruth_trkID",            &fprimaryTruth_trkID,           "primaryTruth_trkID/I");
  fNuBeam->Branch("primaryShowerBestPlane",        &fprimaryShowerBestPlane,       "primaryShowerBestPlane/I");
  fNuBeam->Branch("primaryShowerCharge",           &fprimaryShowerCharge,          "primaryShowerCharge/D");
  fNuBeam->Branch("primaryShowerEnergy",           &fprimaryShowerEnergy,          "primaryShowerEnergy/D");
  fNuBeam->Branch("primaryShowerMIPEnergy",        &fprimaryShowerMIPEnergy,       "primaryShowerMIPEnergy/D");

  fNuBeam->Branch("primaryKineticEnergy",          &fprimaryKineticEnergy,         "primaryKineticEnergy[3]/D");
  fNuBeam->Branch("primaryRange",                  &fprimaryRange,                 "primaryRange[3]/D");
  fNuBeam->Branch("primarynCal",                   &fprimarynCal,                  "primarynCal/I");
  fNuBeam->Branch("primarydQdx",                   &fprimarydQdx,                  "primarydQdx[primarynCal]/D");
  fNuBeam->Branch("primary_cal_pos",               &fprimary_cal_pos,              "primary_cal_pos[primarynCal][3]/D");
  fNuBeam->Branch("primary_cal_pitch",             &fprimary_cal_pitch,            "primary_cal_pitch[primarynCal]/D");
  fNuBeam->Branch("primarydEdx",                   &fprimarydEdx,                  "primarydEdx[primarynCal]/D");
  fNuBeam->Branch("primaryResidualRange",          &fprimaryResidualRange,         "primaryResidualRange[primarynCal]/D");
  fNuBeam->Branch("primaryT0",                     &fprimaryT0,                    "primaryT0/D");
  fNuBeam->Branch("primaryNDaughters",             &fprimaryNDaughters,            "primaryNDaughters/I");

  fNuBeam->Branch("primaryDaughterID",               &fprimaryDaughterID,                "primaryDaughterID[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterIstrack",          &fprimaryDaughterIstrack,           "primaryDaughterIstrack[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterIsshower",         &fprimaryDaughterIsshower,          "primaryDaughterIsshower[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterBDTScore",         &fprimaryDaughterBDTScore,          "primaryDaughterBDTScore[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterCNNScore",         &fprimaryDaughterCNNScore,          "primaryDaughterCNNScore[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterNCollHits",        &fprimaryDaughterNCollHits,         "primaryDaughterNCollHits[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterNHits",            &fprimaryDaughterNHits,             "primaryDaughterNHits[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterEnergy",           &fprimaryDaughterEnergy,            "primaryDaughterEnergy[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterEnergyFromHits",   &fprimaryDaughterEnergyFromHits,    "primaryDaughterEnergyFromHits[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterMomentum",         &fprimaryDaughterMomentum,          "primaryDaughterMomentum[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterEndMomentum",      &fprimaryDaughterEndMomentum,       "primaryDaughterEndMomentum[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterLength",           &fprimaryDaughterLength,            "primaryDaughterLength[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterStartPosition",    &fprimaryDaughterStartPosition,     "primaryDaughterStartPosition[primaryNDaughters][3]/D");
  fNuBeam->Branch("primaryDaughterEndPosition",      &fprimaryDaughterEndPosition,       "primaryDaughterEndPosition[primaryNDaughters][3]/D");
  fNuBeam->Branch("primaryDaughterStartDirection",   &fprimaryDaughterStartDirection,    "primaryDaughterStartDirection[primaryNDaughters][3]/D");
  fNuBeam->Branch("primaryDaughterEndDirection",     &fprimaryDaughterEndDirection,      "primaryDaughterEndDirection[primaryNDaughters][3]/D");
  fNuBeam->Branch("primaryDaughterParentPdg",        &fprimaryDaughterParentPdg,         "primaryDaughterParentPdg[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterGrandparentPdg",   &fprimaryDaughterGrandparentPdg,    "primaryDaughterGrandparentPdg[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterParentID",         &fprimaryDaughterParentID,          "primaryDaughterParentID[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterGrandparentID",    &fprimaryDaughterGrandparentID,     "primaryDaughterGrandparentID[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterParentE",          &fprimaryDaughterParentE,           "primaryDaughterParentE[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterGrandparentE",     &fprimaryDaughterGrandparentE,      "primaryDaughterGrandparentE[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterParentStart",      &fprimaryDaughterParentStart,       "primaryDaughterParentStart[primaryNDaughters][3]/D");
  fNuBeam->Branch("primaryDaughterGrandparentStart", &fprimaryDaughterGrandparentStart,  "primaryDaughterGrandparentStart[primaryNDaughters][3]/D");
  fNuBeam->Branch("primaryDaughterParentEnd",        &fprimaryDaughterParentEnd,         "primaryDaughterParentEnd[primaryNDaughters][3]/D");
  fNuBeam->Branch("primaryDaughterGrandparentEnd",   &fprimaryDaughterGrandparentEnd,    "primaryDaughterGrandparentEnd[primaryNDaughters][3]/D");
  fNuBeam->Branch("primaryDaughterParentEndProcess", &fprimaryDaughterParentEndProcess,  "primaryDaughterParentEndProcess[primaryNDaughters]/I");
  fNuBeam->Branch("primaryDaughterCompleteness",     &fprimaryDaughterCompleteness,      "primaryDaughterCompleteness[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterPurity",           &fprimaryDaughterPurity,            "primaryDaughterPurity[primaryNDaughters]/D");
  fNuBeam->Branch("primaryDaughterHitPosition",      &fprimaryDaughterHitPosition,       ("primaryDaughterHitPosition[primaryNDaughters][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fNuBeam->Branch("primaryDaughterCaloHitPosition",  &fprimaryDaughterCaloHitPosition,   ("primaryDaughterCaloHitPosition[primaryNDaughters][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  fNuBeam->Branch("primaryDaughterHitCharge",        &fprimaryDaughterHitCharge,         ("primaryDaughterHitCharge[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("primaryDaughterHitdEdxArea",      &fprimaryDaughterHitdEdxArea,       ("primaryDaughterHitdEdxArea[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("primaryDaughterHitdQdx",          &fprimaryDaughterHitdQdx,           ("primaryDaughterHitdQdx[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  fNuBeam->Branch("primaryDaughterCaloHitdQdx",      &fprimaryDaughterCaloHitdQdx,       ("primaryDaughterCaloHitdQdx[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  fNuBeam->Branch("primaryDaughterCaloHitPitch",     &fprimaryDaughterCaloHitPitch,      ("primaryDaughterCaloHitPitch[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  fNuBeam->Branch("primaryDaughterOwnHitPitch",      &fprimaryDaughterOwnHitPitch,       ("primaryDaughterOwnHitPitch[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  fNuBeam->Branch("primaryDaughterHitWire",          &fprimaryDaughterHitWire,           ("primaryDaughterHitWire[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/I").c_str());
  fNuBeam->Branch("primaryDaughterHitTime",          &fprimaryDaughterHitTime,           ("primaryDaughterHitTime[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());

  // fNuBeam->Branch("ObjID",                &fObjID,                ("ObjID[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fNuBeam->Branch("ObjIsTrack",           &fObjIsTrack,           ("ObjIsTrack[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fNuBeam->Branch("ObjIsShower",          &fObjIsShower,          ("ObjIsShower[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fNuBeam->Branch("ObjIsPrimaryDaughter", &fObjIsPrimaryDaughter, ("ObjIsPrimaryDaughter[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fNuBeam->Branch("ObjBDTscore",          &fObjBDTscore,          ("ObjBDTscore[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fNuBeam->Branch("ObjCNNscore",          &fObjCNNscore,          ("ObjCNNscore[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fNuBeam->Branch("ObjMomentum",          &fObjMomentum,          ("ObjMomentum[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fNuBeam->Branch("ObjEndMomentum",       &fObjEndMomentum,       ("ObjEndMomentum[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fNuBeam->Branch("ObjLength",            &fObjLength,            ("ObjLength[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fNuBeam->Branch("ObjStartPosition",     &fObjStartPosition,     ("ObjStartPosition[" + std::to_string(NMAXOBJECTS) + "][3]/D").c_str());
  // fNuBeam->Branch("ObjEndPosition",       &fObjEndPosition,       ("ObjEndPosition[" + std::to_string(NMAXOBJECTS) + "][3]/D").c_str());
  // fNuBeam->Branch("ObjStartDirection",    &fObjStartDirection,    ("ObjStartDirection[" + std::to_string(NMAXOBJECTS) + "][3]/D").c_str());
  // fNuBeam->Branch("ObjEndDirection",      &fObjEndDirection,      ("ObjEndDirection[" + std::to_string(NMAXOBJECTS) + "][3]/D").c_str());
  // fNuBeam->Branch("ObjParentPdg",         &fObjMCParentPdg,         ("ObjParentPdg[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fNuBeam->Branch("ObjGrandparentPdg",    &fObjMCGrandparentPdg,    ("ObjGrandparentPdg[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fNuBeam->Branch("ObjParentID",          &fObjMCParentID,          ("ObjParentID[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fNuBeam->Branch("ObjGrandparentID",     &fObjMCGrandparentID,     ("ObjGrandparentID[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // // Object hit information.
  // fNuBeam->Branch("ObjNumHits",           &fObjNumHits,           ("ObjNumHits[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fNuBeam->Branch("ObjHitPosition",       &fObjHitPosition,       ("ObjHitPosition[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fNuBeam->Branch("ObjHitTime",           &fObjHitTime,           ("ObjHitTime[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("ObjHitWire",           &fObjHitWire,           ("ObjHitWire[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/I").c_str());
  // fNuBeam->Branch("ObjHitCharge",         &fObjHitCharge,         ("ObjHitCharge[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("ObjHitPitch",          &fObjHitPitch,          ("ObjHitPitch[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("ObjHitdEdx",           &fObjHitdEdx,           ("ObjHitdEdx[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("ObjHitdQdx",           &fObjHitdQdx,           ("ObjHitdQdx[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());

  // MC pi0 variables
  fNuBeam->Branch("NMCPi0s",                       &fNMCPi0s,                      "NMCPi0s/I");
  fNuBeam->Branch("MCPi0ID",                       &fMCPi0ID,                      "MCPi0ID[NMCPi0s]/I");
  fNuBeam->Branch("MCPi0Energy",                   &fMCPi0Energy,                  "MCPi0Energy[NMCPi0s]/D");
  fNuBeam->Branch("MCPi0StartPosition",            &fMCPi0StartPosition,           "MCPi0StartPosition[NMCPi0s][3]/D");
  fNuBeam->Branch("MCPi0EndPosition",              &fMCPi0EndPosition,             "MCPi0EndPosition[NMCPi0s][3]/D");
  fNuBeam->Branch("MCPhoton1ID",                   &fMCPhoton1ID,                  "MCPhoton1ID[NMCPi0s]/I");
  fNuBeam->Branch("MCPhoton2ID",                   &fMCPhoton2ID,                  "MCPhoton2ID[NMCPi0s]/I");
  fNuBeam->Branch("MCPhoton1Energy",               &fMCPhoton1Energy,              "MCPhoton1Energy[NMCPi0s]/D");
  fNuBeam->Branch("MCPhoton2Energy",               &fMCPhoton2Energy,              "MCPhoton2Energy[NMCPi0s]/D");
  fNuBeam->Branch("MCPhoton1StartPosition",        &fMCPhoton1StartPosition,       "MCPhoton1StartPosition[NMCPi0s][3]/D");
  fNuBeam->Branch("MCPhoton2StartPosition",        &fMCPhoton2StartPosition,       "MCPhoton2StartPosition[NMCPi0s][3]/D");
  fNuBeam->Branch("MCPhoton1EndPosition",          &fMCPhoton1EndPosition,         "MCPhoton1EndPosition[NMCPi0s][3]/D");
  fNuBeam->Branch("MCPhoton2EndPosition",          &fMCPhoton2EndPosition,         "MCPhoton2EndPosition[NMCPi0s][3]/D");
  // // Reco hits related to photons
  // fNuBeam->Branch("MCPhoton1NumHits",              &fMCPhoton1NumHits,             ("MCPhoton1NumHits[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("MCPhoton2NumHits",              &fMCPhoton2NumHits,             ("MCPhoton2NumHits[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("MCPhoton1_hit_w",               &fMCPhoton1_hit_w,              ("MCPhoton1_hit_w[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/I").c_str());
  // fNuBeam->Branch("MCPhoton2_hit_w",               &fMCPhoton2_hit_w,              ("MCPhoton2_hit_w[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/I").c_str());
  // fNuBeam->Branch("MCPhoton1_hit_q",               &fMCPhoton1_hit_q,              ("MCPhoton1_hit_q[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("MCPhoton2_hit_q",               &fMCPhoton2_hit_q,              ("MCPhoton2_hit_q[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("MCPhoton1_hit_t",               &fMCPhoton1_hit_t,              ("MCPhoton1_hit_t[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("MCPhoton2_hit_t",               &fMCPhoton2_hit_t,              ("MCPhoton2_hit_t[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("MCPhoton1_hit_pos",             &fMCPhoton1_hit_pos,            ("MCPhoton1_hit_pos[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fNuBeam->Branch("MCPhoton2_hit_pos",             &fMCPhoton2_hit_pos,            ("MCPhoton2_hit_pos[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());

  // // Shower variables
  // fNuBeam->Branch("Shower1ID",                     &fShower1ID,                    ("Shower1ID[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("Shower2ID",                     &fShower2ID,                    ("Shower2ID[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("Shower1StartPosition",          &fShower1StartPosition,         ("Shower1StartPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fNuBeam->Branch("Shower2StartPosition",          &fShower2StartPosition,         ("Shower2StartPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fNuBeam->Branch("Shower1Direction",              &fShower1Direction,             ("Shower1Direction[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fNuBeam->Branch("Shower2Direction",              &fShower2Direction,             ("Shower2Direction[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fNuBeam->Branch("Shower1Length",                 &fShower1Length,                ("Shower1Length[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2Length",                 &fShower2Length,                ("Shower2Length[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower1Completeness",           &fShower1Completeness,          ("Shower1Completeness[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2Completeness",           &fShower2Completeness,          ("Shower2Completeness[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower1Purity",                 &fShower1Purity,                ("Shower1Purity[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2Purity",                 &fShower2Purity,                ("Shower2Purity[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // // Reco hits related to showers
  // fNuBeam->Branch("Shower1NumHits",                &fShower1NumHits,               ("Shower1NumHits[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("Shower2NumHits",                &fShower2NumHits,               ("Shower2NumHits[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("Shower1_cnn_sc",                &fShower1_cnn_sc,               ("Shower1_cnn_sc[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2_cnn_sc",                &fShower2_cnn_sc,               ("Shower2_cnn_sc[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("Shower1_cal_pos",               &fShower1_cal_pos,              ("Shower1_cal_pos[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fNuBeam->Branch("Shower2_cal_pos",               &fShower2_cal_pos,              ("Shower2_cal_pos[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fNuBeam->Branch("Shower1_cal_E",                 &fShower1_cal_E,                ("Shower1_cal_E[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2_cal_E",                 &fShower2_cal_E,                ("Shower2_cal_E[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower1Energy",                 &fShower1Energy,                ("Shower1Energy[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2Energy",                 &fShower2Energy,                ("Shower2Energy[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower1EnergyFromHits",         &fShower1EnergyFromHits,        ("Shower1EnergyFromHits[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2EnergyFromHits",         &fShower2EnergyFromHits,        ("Shower2EnergyFromHits[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fNuBeam->Branch("Shower1HasBeamParent",          &fShower1HasBeamParent,         ("Shower1HasBeamParent[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("Shower2HasBeamParent",          &fShower2HasBeamParent,         ("Shower2HasBeamParent[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("Shower1HasBeamGrandparent",     &fShower1HasBeamGrandparent,    ("Shower1HasBeamGrandparent[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("Shower2HasBeamGrandparent",     &fShower2HasBeamGrandparent,    ("Shower2HasBeamGrandparent[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fNuBeam->Branch("Shower1_cal_pitch",             &fShower1_cal_pitch,            ("Shower1_cal_pitch[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2_cal_pitch",             &fShower2_cal_pitch,            ("Shower2_cal_pitch[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("Shower1_cal_dEdx",              &fShower1_cal_dEdx,             ("Shower1_cal_dEdx[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2_cal_dEdx",              &fShower2_cal_dEdx,             ("Shower2_cal_dEdx[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("Shower1_cal_dQdx",              &fShower1_cal_dQdx,             ("Shower1_cal_dQdx[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fNuBeam->Branch("Shower2_cal_dQdx",              &fShower2_cal_dQdx,             ("Shower2_cal_dQdx[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
}

void protoana::DUNEPizeroAnaTree::analyze(art::Event const & evt){
  // Initialise tree parameters
  Initialise();
  fRun = evt.run();
  fSubRun = evt.subRun();
  fevent  = evt.id().event();
  art::Timestamp ts = evt.time();
  if (ts.timeHigh() == 0){
    TTimeStamp ts2(ts.timeLow());
    fTimeStamp = ts2.AsDouble();
  }
  else{
    TTimeStamp ts2(ts.timeHigh(), ts.timeLow());
    fTimeStamp = ts2.AsDouble();
  }

  // // Get number of active fembs
  // if(!evt.isRealData()){
  //   for(int k=0; k < 6; k++)
  //     fNactivefembs[k] = 20;
  // }
  // else{
  //   for(int k=0; k < 6; k++)
  //     fNactivefembs[k] = dataUtil.GetNActiveFembsForAPA(evt, k);
  // }

  //check for reco pandora stuff
  art::Handle<std::vector<recob::PFParticle>> recoParticleHandle;
  if( !evt.getByLabel(fPFParticleTag,recoParticleHandle) ) return;

  // Get all of the PFParticles, by default from the "pandora" product
  // auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  for(const recob::PFParticle* particle : pfParticles){

    FillPrimaryPFParticle(evt, particle);
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackTag);
    // std::cout << "Primary PFParticle ID: " << particle->Self() << '\n';

    fprimaryVertex[0] = vtx.X(); fprimaryVertex[1] = vtx.Y(); fprimaryVertex[2] = vtx.Z();
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackTag);

    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
    break;
  }

  // // Look at all the reconstructed particles that are fully contained in APA3.
  // FillAPA3Object(evt, recoParticleHandle);

  bool beamTriggerEvent = false;
  if(!evt.isRealData()){
    // // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    // auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    // auto genTruth = (*mcTruths)[0];
    // Define and fill a handle to point to a vector of the MCParticles
    art::Handle<std::vector<simb::MCParticle>> MCParticleHandle;
    if (!evt.getByLabel(fSimulationTag, MCParticleHandle)) {
      // Handle no simb::MCParticles.
      throw cet::exception("DUNEPizeroAnaTree")
          << " No simb::MCParticle objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    std::map<int, const simb::MCParticle*> MCPmap;
    for(const simb::MCParticle& mcp : *MCParticleHandle) {
      MCPmap[mcp.TrackId()] = &mcp;
    }

    // Also get the reconstructed beam information in the MC - TO DO
    // const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    // Get first primary MCParticle.
    std::vector<const simb::MCParticle*> primNus;
    const simb::MCParticle* geantGoodParticle = 0x0;
    // for(const simb::MCParticle& mcp : *MCParticleHandle) {
    //   const int gPDG = abs(mcp.PdgCode());
    //   const bool particleIsNeutrino = gPDG == 12 || gPDG == 14 || gPDG == 16;
    //   if(mcp.Process() == "primary" && particleIsNeutrino) {
    //     geantGoodParticle = &mcp;
    //     // std::cout << mcp.PdgCode() << '\n';
    //     break;
    //   }
    // }

    // // Get list of the g4 particles. plist should be a std::map< int, simb::MCParticle* >
    art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
    const sim::ParticleList & plist = pi_serv->ParticleList();
    for(auto const part : plist){
      // std::cout << "PDG: " << part.second->PdgCode() << ", process: " << part.second->Process() << '\n';
      if(part.second->Process() == "primary" && (part.second->PdgCode() == 12 || part.second->PdgCode() == 14 || part.second->PdgCode() == 16)){
        primNus.push_back(part.second);
        geantGoodParticle = part.second;
        // std::cout << geantGoodParticle << '\n';
        break;
      }
    }

    // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    if(geantGoodParticle != 0x0) {
      beamTriggerEvent = true;
      fbeamtrigger       = 12;
      fbeamtrackPos[0]   = geantGoodParticle->Vx();
      fbeamtrackPos[1]   = geantGoodParticle->Vy();
      fbeamtrackPos[2]   = geantGoodParticle->Vz();
      fbeamtrackEndPos[0]= geantGoodParticle->EndX();
      fbeamtrackEndPos[1]= geantGoodParticle->EndY();
      fbeamtrackEndPos[2]= geantGoodParticle->EndZ();
      fbeamtrackMomentum = geantGoodParticle->P();
      fbeamtrackEnergy   = geantGoodParticle->E();
      fbeamtrackPdg      = geantGoodParticle->PdgCode();
      fbeamtrackTime     = geantGoodParticle->T();
      fbeamtrackID       = geantGoodParticle->TrackId();
      fbeamtrackEndProcess.push_back(geantGoodParticle->EndProcess());
      fbeam_ntrjPoints = std::min(geantGoodParticle->NumberTrajectoryPoints(), (unsigned)1000);
      fbeamtrackNDaughters = geantGoodParticle->NumberDaughters();
      for( size_t i=0; i<geantGoodParticle->NumberTrajectoryPoints() && i<1000; ++i){
         fbeamtrackPos_at[i][0] = geantGoodParticle->Position(i).X();
         fbeamtrackPos_at[i][1] = geantGoodParticle->Position(i).Y();
         fbeamtrackPos_at[i][2] = geantGoodParticle->Position(i).Z();
         fbeamtrackMom_at[i][0] = geantGoodParticle->Momentum(i).Px();
         fbeamtrackMom_at[i][1] = geantGoodParticle->Momentum(i).Py();
         fbeamtrackMom_at[i][2] = geantGoodParticle->Momentum(i).Pz();
         fbeamtrackMom_at[i][3] = geantGoodParticle->Momentum(i).E();
      }

      std::cout << "Primary beam particle PDG code: " << fbeamtrackPdg << " and E: "
                << fbeamtrackEnergy << '\n';
      // std::cout << fbeamtrackEndPos[0] << ',' << fbeamtrackEndPos[1] << ',' << fbeamtrackEndPos[2] << "\n\n";

      // art::Handle<std::vector<recob::Vertex>> VertexHandle;
      // if (!evt.getByLabel(fPFParticleTag, VertexHandle)) {
      //   // Handle no recob::Vertexs.
      //   throw cet::exception("DUNEPizeroAnaTree")
      //       << " No recob::Vertex objects in this event - "
      //       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      // }
      // for(const recob::Vertex& vtx : *VertexHandle) {
      //   std::cout << vtx.position().x() << ',' << vtx.position().y() << ',' << vtx.position().z() << '\n';
      // }

      std::vector<const simb::MCParticle*> pi0s;

      // Search for daughter pi0s of the primary particle.
      if(geantGoodParticle->NumberDaughters() > 0) {
        std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
        std::cout << "Primary has " << geantGoodParticle->NumberDaughters() << " daughters.\n";
      }

      std::cout << "ProcessMap size: " << geantGoodParticle->Trajectory().TrajectoryProcesses().size() << '\n';
      for(const std::pair<size_t, unsigned char>& p : geantGoodParticle->Trajectory().TrajectoryProcesses()) {
        std::cout << p.first << "/" << geantGoodParticle->NumberTrajectoryPoints() << " , "
                  << geantGoodParticle->Trajectory().KeyToProcess(p.second) << '\n';
      }

      std::cout << "End process: " << geantGoodParticle->EndProcess() << ".\n";

      // Find distance between neutrino end and close other particles.
      const TLorentzVector nuEnd = geantGoodParticle->EndPosition();
      for(auto const part : plist) {
        const simb::MCParticle& mcp = *(part.second);
        const double dist = (mcp.Position() - nuEnd).Vect().Mag();
        if(dist < 100) {
          std::cout << "MCP ID " << mcp.TrackId() << ", PDG " << mcp.PdgCode()
                    << ", distance " << dist << '\n';
        }
      }

      // std::cout << "Primary has pi0 daughters:\n";
      // std::vector<const simb::MCParticle*> daughters;
      // daughters.push_back(geantGoodParticle);
      // while(daughters.size()) {
      //   const simb::MCParticle* parent = daughters[0];
      //   for(int di = 0; di < parent->NumberDaughters(); ++di) {
      //     const int daughter_ID = parent->Daughter(di);
      //     const simb::MCParticle* daughter = MCPmap[daughter_ID];
      //     daughters.push_back(daughter);
      //     if(daughter->PdgCode() == 111) {
      //       std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
      //       std::cout << "  ID " << daughter_ID << ", PDG " << daughter->PdgCode()
      //                 << ", P " << daughter->P() << '\n';
      //       pi0s.push_back(daughter);
      //     }
      //     std::cout << daughter->PdgCode() << ',' << daughter->E() << '\n';
      //   }
      //   daughters.erase(daughters.begin());
      // }

      // Record MC and shower information
      setPiZeroInfo(evt, pi0s);

    } //geantGoodParticle
  } // MC
  else { // data
    std::cout << "Is this the future? We don't have DUNE FD data yet in 2020.\n";
  }// data

  // Fill tree if event came from beam trigger.
  if(beamTriggerEvent)   fNuBeam->Fill();

  // // Put some space between events.
  // std::cout << "\n\n";
}

// -----------------------------------------------------------------------------
void protoana::DUNEPizeroAnaTree::endJob(){

}

// -----------------------------------------------------------------------------
void protoana::DUNEPizeroAnaTree::FillPrimaryPFParticle(art::Event const & evt, const recob::PFParticle* particle){

  // Pandora's BDT beam-cosmic score
  fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);

  // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
  // of this particle might be more helpful. These return null pointers if not track-like / shower-like
  const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle, evt,fPFParticleTag,fTrackTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track>>(fTrackTag);
  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower>>(fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromTracks(recoTracks,evt,fTrackTag);
  // Set some common variables that have different access methods.
  std::vector<art::Ptr<recob::Hit>> pfpHits;
  // std::vector<anab::Calorimetry> calovector;
  const simb::MCParticle* mcparticle = NULL;
  if(thisTrack != 0x0) {
    pfpHits = findHitsFromTracks.at(thisTrack->ID());
    // calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackTag, fCalorimetryTag);
    mcparticle = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackTag);
  } else if(thisShower != 0x0) {
    pfpHits = findHitsFromShowers.at(thisShower->ID());
    // calovector = showerUtil.GetRecoShowerCalorimetry(*thisShower, evt, fShowerTag, fShowerCaloTag);
    mcparticle = truthUtil.GetMCParticleFromRecoShower(*thisShower, evt, fShowerTag);
  }
  // NHits associated with this pfParticle.
  fprimaryNHits = pfpHits.size();

  if(mcparticle){
    fprimaryTruth_pdg = mcparticle->PdgCode();
    fprimaryTruth_trkID = mcparticle->TrackId();
    fprimaryTruth_vtx[0] = mcparticle->Vx();
    fprimaryTruth_vtx[1] = mcparticle->Vy();
    fprimaryTruth_vtx[2] = mcparticle->Vx();
    fprimaryTruth_E = mcparticle->E();
    if( fbeamtrackID != -999 && fbeamtrackID == fprimaryTruth_trkID ) {
      fprimaryIsBeamparticle = 1;
    }
  }

  // // Calorimetry only collection plane
  // for(size_t k = 0; k < calovector.size() && k<3; k++){
  //   int plane = calovector[k].PlaneID().Plane;
  //   if(plane !=2 ) continue;
  //   fprimaryKineticEnergy[plane]      = calovector[k].KineticEnergy();
  //   fprimaryRange[plane]              = calovector[k].Range();
  //   fprimarynCal = std::min(calovector[k].dEdx().size(), (size_t)NMAXHITS);
  //   for(size_t l=0; l<calovector[k].dEdx().size() && l<NMAXHITS; ++l){
  //     fprimarydEdx[l]= calovector[k].dEdx()[l];
  //     fprimarydQdx[l]= calovector[k].dQdx()[l];
  //     fprimaryResidualRange[l]= calovector[k].ResidualRange()[l];
  //     const auto &pos=(calovector[k].XYZ())[l];
  //     fprimary_cal_pos[l][0] = pos.X();
  //     fprimary_cal_pos[l][1] = pos.Y();
  //     fprimary_cal_pos[l][2] = pos.Z();
  //     fprimary_cal_pitch[l] =calovector[k].TrkPitchVec()[l];
  //   }
  // }

  // Get the CNN score for the primary particle.
  anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );
  fprimaryCNNScore = 0;
  for(const art::Ptr<recob::Hit>& hit : pfpHits) {
    std::array<float,4> cnn_out = hitResults.getOutput( hit );
    const double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ];
    const double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh;
    fprimaryCNNScore += cnn_score;
  }
  fprimaryCNNScore /= fprimaryNHits;
  // std::cout << "Got primary CNN score " << fprimaryCNNScore << ", true PDG code: " << fprimaryTruth_pdg << '\n';
  // std::cout << "Pandora track: " << thisTrack << ", shower: " << thisShower << '\n';

  // Get the T0 for this pfParticle
  std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
  if(!pfT0vec.empty())
    fprimaryT0 = pfT0vec[0].Time();

  if(thisTrack != 0x0){

    fprimaryIstrack                    = 1;
    fprimaryIsshower                   = 0;
    fprimaryID                         = thisTrack->ParticleId();
    fprimaryTheta                      = thisTrack->Theta();
    fprimaryPhi                        = thisTrack->Phi();
    fprimaryLength                     = thisTrack->Length();
    fprimaryMomentum                   = thisTrack->StartMomentum();
    fprimaryEndMomentum                = thisTrack->EndMomentum();
    fprimaryEndPosition[0]             = thisTrack->Trajectory().End().X();
    fprimaryEndPosition[1]             = thisTrack->Trajectory().End().Y();
    fprimaryEndPosition[2]             = thisTrack->Trajectory().End().Z();
    fprimaryStartPosition[0]           = thisTrack->Trajectory().Start().X();
    fprimaryStartPosition[1]           = thisTrack->Trajectory().Start().Y();
    fprimaryStartPosition[2]           = thisTrack->Trajectory().Start().Z();
    fprimaryEndDirection[0]            = thisTrack->Trajectory().EndDirection().X();
    fprimaryEndDirection[1]            = thisTrack->Trajectory().EndDirection().Y();
    fprimaryEndDirection[2]            = thisTrack->Trajectory().EndDirection().Z();
    fprimaryStartDirection[0]          = thisTrack->Trajectory().StartDirection().X();
    fprimaryStartDirection[1]          = thisTrack->Trajectory().StartDirection().Y();
    fprimaryStartDirection[2]          = thisTrack->Trajectory().StartDirection().Z();

  } // end is track
  else if(thisShower != 0x0){

    fprimaryIstrack                     = 0;
    fprimaryIsshower                    = 1;
    fprimaryID                          = thisShower->ID();
    fprimaryLength                      = thisShower->Length();
    fprimaryShowerBestPlane             = thisShower->best_plane();
    fprimaryOpeningAngle                = thisShower->OpenAngle();
    fprimaryStartPosition[0]            = thisShower->ShowerStart().X();
    fprimaryStartPosition[1]            = thisShower->ShowerStart().Y();
    fprimaryStartPosition[2]            = thisShower->ShowerStart().Z();
    fprimaryStartDirection[0]           = thisShower->Direction().X();
    fprimaryStartDirection[1]           = thisShower->Direction().Y();
    fprimaryStartDirection[2]           = thisShower->Direction().Z();

  } // end is shower

  // Record children of primary track
  art::Handle<std::vector<recob::PFParticle>> recoParticleHandle;
  if( !evt.getByLabel(fPFParticleTag,recoParticleHandle) ) return;

  fprimaryNDaughters = std::min(particle->NumDaughters(), NMAXDAUGHTERS);
  for(int di = 0; di < fprimaryNDaughters; ++di) {
    fprimaryDaughterID[di] = particle->Daughter(di);
    const recob::PFParticle* daughter = &recoParticleHandle->at(fprimaryDaughterID[di]);
    if(daughter == 0x0) continue;

    const recob::Track* daughterTrack   = pfpUtil.GetPFParticleTrack(*daughter, evt,fPFParticleTag,fTrackTag);
    const recob::Shower* daughterShower = pfpUtil.GetPFParticleShower(*daughter,evt,fPFParticleTag,fShowerTag);

    // Pandora's BDT beam-cosmic score
    fprimaryDaughterBDTScore[di] = (double)pfpUtil.GetBeamCosmicScore(*daughter,evt,fPFParticleTag);
    // Hits and calorimetry associated with this pfParticle
    std::vector<art::Ptr<recob::Hit>> pfpDHits;
    // std::vector<anab::Calorimetry> dcalovector;
    if(daughterTrack != 0x0) {
      pfpDHits = findHitsFromTracks.at( daughterTrack->ID() );
      // dcalovector = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackTag, fCalorimetryTag);
    } else if(daughterShower != 0x0) {
      pfpDHits = findHitsFromShowers.at( daughterShower->ID() );
      // dcalovector = showerUtil.GetRecoShowerCalorimetry(*daughterShower, evt, fShowerTag, fShowerCaloTag);
    }
    // if(dcalovector.size() != 3 && fVerbose > 0) {
    //   std::cerr << "WARNING::Calorimetry vector size for track is = " << dcalovector.size() << std::endl;
    // }

    // // Loop over calorimetry.
    // std::map<int, double> calo_hit_dQdx;
    // std::map<int, double> calo_hit_pitch;
    // for(size_t k = 0; k < dcalovector.size() && k<3; k++){
    //   int plane = dcalovector[k].PlaneID().Plane;
    //   if(plane !=2 ) continue;
    //   for(size_t l=0; l<dcalovector[k].dEdx().size() && l<NMAXHITS; ++l){
    //     // fprimaryDaughterHitdEdx[di][l]= dcalovector[k].dEdx()[l];
    //     // fprimaryDaughterCaloHitdQdx[di][l]= dcalovector[k].dQdx()[l];
    //     calo_hit_dQdx[dcalovector[k].TpIndices()[l]] = dcalovector[k].dQdx()[l];
    //     calo_hit_pitch[dcalovector[k].TpIndices()[l]] = dcalovector[k].TrkPitchVec()[l];
    //     // tmp_dQdx.push_back(fprimaryDaughterCaloHitdQdx[di][l]);
    //     // fprimaryDaughterResidualRange[l]= dcalovector[k].ResidualRange()[l];
    //     const geo::Point_t &pos=(dcalovector[k].XYZ())[l];
    //     if(daughterShower != 0) {
    //       const TVector3 hitvec(pos.X(), pos.Y(), pos.Z());
    //     }
    //     // fprimaryDaughterCaloHitPosition[di][l][0] = pos.X();
    //     // fprimaryDaughterCaloHitPosition[di][l][1] = pos.Y();
    //     // fprimaryDaughterCaloHitPosition[di][l][2] = pos.Z();
    //     // fprimaryDaughterHitPitch[di][l] = dcalovector[k].TrkPitchVec()[l];
    //   } // Loop over calorimetry hits.
    // } // Loop over calorimetry.

    // Loop over hits, record Aidan's CNN shower-like score.
    art::FindManyP<recob::SpacePoint> spFromHits(pfpDHits, evt, fHitTag);
    fprimaryDaughterCNNScore[di] = 0;
    unsigned rec_hiti = 0;
    for( unsigned i = 0; i < pfpDHits.size() && rec_hiti < NMAXHITS; ++i){
      if(pfpDHits[i]->WireID().Plane != 2) continue;

      const geo::WireGeo* pwire = fGeometryService->WirePtr(pfpDHits[i]->WireID());
      TVector3 xyzWire = pwire->GetCenter<TVector3>();
      fprimaryDaughterHitWire[di][rec_hiti] = pfpDHits[i]->WireID().Wire;
      fprimaryDaughterHitTime[di][rec_hiti] = pfpDHits[i]->PeakTime();
      fprimaryDaughterHitPosition[di][rec_hiti][0] = detprop->ConvertTicksToX(pfpDHits[i]->PeakTime(),pfpDHits[i]->WireID().Plane,pfpDHits[i]->WireID().TPC,0);
      fprimaryDaughterHitPosition[di][rec_hiti][2] = xyzWire.Z();
      std::vector<art::Ptr<recob::SpacePoint>> sps = spFromHits.at(i);
      if(!sps.empty()) {
        fprimaryDaughterHitPosition[di][rec_hiti][1] = sps[0]->XYZ()[1];
      }
      TVector3 hitpos(fprimaryDaughterHitPosition[di][rec_hiti][0],
                      fprimaryDaughterHitPosition[di][rec_hiti][1],
                      fprimaryDaughterHitPosition[di][rec_hiti][2]);

      std::array<float,4> cnn_out = hitResults.getOutput( pfpDHits[i] );
      const double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ];
      const double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh;
      fprimaryDaughterCNNScore[di] += cnn_score;
      // Record general hit information.
      fprimaryDaughterHitCharge[di][rec_hiti] = pfpDHits[i]->Integral();
      // Determine own hit pitch with SCE corrections.
      const double wirePitch = fGeometryService->WirePitch(pfpDHits[i]->WireID().planeID());
      TVector3 objdir(0,0,0);
      if(daughterTrack != 0x0) {
        objdir = {daughterTrack->Trajectory().StartDirection().X(),
                  daughterTrack->Trajectory().StartDirection().Y(),
                  daughterTrack->Trajectory().StartDirection().Z()};
      } else if(daughterShower != 0x0) {
        objdir = daughterShower->Direction();
      }
      TVector3 hitdir = objdir.Unit();
      if(objdir.Mag2() != 0 && hitpos[1] != -999.0) {
        // Here we have the object's direction and full hit position.
        // Set hit end to point on next wire in object's direction.
        objdir = objdir.Unit();
        TVector3 hitend = hitpos + objdir * (wirePitch/objdir.Z());
        if(fSCEService->EnableCalSpatialSCE()) {
          // Do SCE corrections on hitpos and hitend to get proper hit direction.
          // std::cout << "Doing SCE corrections.\n";
          const geo::Vector_t poscorr = fSCEService->GetCalPosOffsets(geo::Point_t(hitpos), pfpDHits[i]->WireID().TPC);
          const geo::Vector_t endcorr = fSCEService->GetCalPosOffsets(geo::Point_t(hitend), pfpDHits[i]->WireID().TPC);
          hitpos += TVector3(poscorr.X(), poscorr.Y(), poscorr.Z());
          hitend += TVector3(endcorr.X(), endcorr.Y(), endcorr.Z());
        }
        hitdir = (hitend-hitpos).Unit();
      } //else if(objdir.Mag2() != 0) {
      //   std::cout << "Hit position incomplete: " << hitpos[0] << ',' << hitpos[1] << ',' << hitpos[2] << '\n';
      // } else if(hitpos[1] != -999) {
      //   std::cout << "No object direction.\n";
      // } else { std::cout << "??????\n"; }
      fprimaryDaughterOwnHitPitch[di][rec_hiti] = wirePitch / abs(hitdir.Z());
      // fprimaryDaughterCaloHitdQdx[di][rec_hiti] = calo_hit_dQdx[pfpDHits[i].key()];
      // fprimaryDaughterCaloHitPitch[di][rec_hiti] = calo_hit_pitch[pfpDHits[i].key()];
      // std::cout << "Bare dx: " << wirePitch / abs(objdir.Z())
      //           << ", own dx: " << fprimaryDaughterOwnHitPitch[di][rec_hiti]
      //           << ", calo dx: " << fprimaryDaughterCaloHitPitch[di][rec_hiti] << '\n';
      ++rec_hiti;
    } // Loop over hits.
    fprimaryDaughterNCollHits[di] = rec_hiti;
    fprimaryDaughterNHits[di] = pfpDHits.size();
    fprimaryDaughterCNNScore[di] /= rec_hiti;

    if(daughterTrack != 0x0) {
      fprimaryDaughterIstrack[di] = 1;
      fprimaryDaughterIsshower[di] = 0;

      fprimaryDaughterMomentum[di]           = daughterTrack->StartMomentum();
      fprimaryDaughterEndMomentum[di]        = daughterTrack->EndMomentum();
      fprimaryDaughterLength[di]             = daughterTrack->Length();
      fprimaryDaughterStartPosition[di][0]   = daughterTrack->Trajectory().Start().X();
      fprimaryDaughterStartPosition[di][1]   = daughterTrack->Trajectory().Start().Y();
      fprimaryDaughterStartPosition[di][2]   = daughterTrack->Trajectory().Start().Z();
      fprimaryDaughterEndPosition[di][0]     = daughterTrack->Trajectory().End().X();
      fprimaryDaughterEndPosition[di][1]     = daughterTrack->Trajectory().End().Y();
      fprimaryDaughterEndPosition[di][2]     = daughterTrack->Trajectory().End().Z();
      fprimaryDaughterStartDirection[di][0]  = daughterTrack->Trajectory().StartDirection().X();
      fprimaryDaughterStartDirection[di][1]  = daughterTrack->Trajectory().StartDirection().Y();
      fprimaryDaughterStartDirection[di][2]  = daughterTrack->Trajectory().StartDirection().Z();
      fprimaryDaughterEndDirection[di][0]    = daughterTrack->Trajectory().EndDirection().X();
      fprimaryDaughterEndDirection[di][1]    = daughterTrack->Trajectory().EndDirection().Y();
      fprimaryDaughterEndDirection[di][2]    = daughterTrack->Trajectory().EndDirection().Z();
    } else if(daughterShower != 0x0) {
      fprimaryDaughterIstrack[di] = 0;
      fprimaryDaughterIsshower[di] = 1;

      std::vector<const recob::Hit*> shHits = showerUtil.GetRecoShowerHits(*daughterShower, evt, fShowerTag);
      fprimaryDaughterEnergyFromHits[di]     = showerUtil.EstimateEnergyFromHitCharge(shHits, fCalorimetryAlg)[2];
      fprimaryDaughterLength[di]             = daughterShower->Length();
      fprimaryDaughterStartDirection[di][0]   = daughterShower->Direction().X();
      fprimaryDaughterStartDirection[di][1]   = daughterShower->Direction().Y();
      fprimaryDaughterStartDirection[di][2]   = daughterShower->Direction().Z();
      fprimaryDaughterStartPosition[di][0]    = daughterShower->ShowerStart().X();
      fprimaryDaughterStartPosition[di][1]    = daughterShower->ShowerStart().Y();
      fprimaryDaughterStartPosition[di][2]    = daughterShower->ShowerStart().Z();
    }

    if(!evt.isRealData()) {
      // Attempt to get the (grand)parent truth information of this daughter.
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      art::ServiceHandle<cheat::BackTrackerService> bt_serv;
      std::unordered_map<int, double> hitMCPEmap;
      for(const art::Ptr<recob::Hit>& hitp : pfpDHits) {
        for(const sim::TrackIDE& ide : bt_serv->HitToTrackIDEs(hitp)) {
          hitMCPEmap[abs(ide.trackID)] += ide.energy;
        }
      }
      int parentID = -1;
      double parentE = -1;
      for(const std::pair<int, double> p : hitMCPEmap) {
        if(p.second > parentE) {
          parentID = p.first;
          parentE = p.second;
        }
      }

      // const simb::MCParticle* parent = truthUtil.GetMCParticleFromReco(*daughter, evt, fPFParticleTag)
      const simb::MCParticle* parent = pi_serv->TrackIdToParticle_P(parentID);
      // Fill MC information if it was found.
      if(parent != 0x0) {
        fprimaryDaughterParentPdg[di]          = parent->PdgCode();
        // std::cout << "True parent PDG: " << fprimaryDaughterParentPdg[di] << '\n';
        fprimaryDaughterParentID[di]          = parent->TrackId();
        // std::cout << "Daughter parent: " << fprimaryDaughterParentID[di] << '\n';
        fprimaryDaughterParentE[di]          = parent->E();
        // std::cout << "True parent E: " << fprimaryDaughterParentE[di] << '\n';
        if(parent->EndProcess() == "conv") {
          fprimaryDaughterParentEndProcess[di] = 0;
        } else if(parent->EndProcess() == "phot") {
          fprimaryDaughterParentEndProcess[di] = 1;
        }
        parent->Position().Vect().GetXYZ(fprimaryDaughterParentStart[di]);
        parent->EndPosition().Vect().GetXYZ(fprimaryDaughterParentEnd[di]);
        if(parent->Mother() != 0) {
          const simb::MCParticle* gparent = pi_serv->TrackIdToParticle_P(parent->Mother());
          fprimaryDaughterGrandparentPdg[di]   =  gparent->PdgCode();
          fprimaryDaughterGrandparentID[di]   =  gparent->TrackId();
          fprimaryDaughterGrandparentE[di]   =  gparent->E();
          // std::cout << "True grandparent PDG: " << fprimaryDaughterGrandparentPdg[di] << '\n';
          gparent->Position().Vect().GetXYZ(fprimaryDaughterGrandparentStart[di]);
          gparent->EndPosition().Vect().GetXYZ(fprimaryDaughterGrandparentEnd[di]);
        }

        // Record completeness and purity.
        // const unsigned numMChits = truthUtil.GetMCParticleHits(*parent, evt, fHitTag).size();
        // unsigned sharedHits = truthUtil.GetSharedHits(*parent, *daughter, evt, fPFParticleTag, false).size() + truthUtil.GetSharedHits(*parent, *daughter, evt, fPFParticleTag, true).size();
        // if(daughterTrack != 0x0) {
        //   sharedHits = truthUtil.GetSharedHits(*parent, *daughter, evt, fPFParticleTag, false).size();
        // } else if(daughterShower != 0x0) {
        //   sharedHits = truthUtil.GetSharedHits(*parent, *daughter, evt, fPFParticleTag, true).size();
        // }
        // Find number of shared hits.
        unsigned sharedHits = 0;
        for(const art::Ptr<recob::Hit>& hitp : pfpDHits) {
          for(const sim::TrackIDE& ide : bt_serv->HitToTrackIDEs(hitp)) {
            if(abs(ide.trackID) == parentID) {
              ++sharedHits;
              break;
            }
          }
        }
        // Find number of MC hits from all hits in detector.
        art::Handle<std::vector<recob::Hit>> hitHandle;
        if(!evt.getByLabel(fHitTag, hitHandle)) {
          std::cout << "Could not find hits in event.\n";
        }
        unsigned numMChits = 0;
        for(const recob::Hit& hit : *hitHandle) {
          for(const int trackId : bt_serv->HitToTrackIds(hit)) {
            if(abs(trackId) == abs(parentID)) {
              ++numMChits;
              break;
            }
          }
        }
        fprimaryDaughterCompleteness[di] = numMChits > 0? (double)sharedHits/numMChits : -999.0;
        fprimaryDaughterPurity[di] = pfpDHits.size() > 0? (double)sharedHits/pfpDHits.size() : -999.0;
        // std::cout << "SHARED HITS: " << sharedHits << ", PDG CODE: " << daughter->PdgCode() << '\n';
        // std::cout << "COMP: " << fprimaryDaughterCompleteness[di] << ", PUR: " << fprimaryDaughterPurity[di] << '\n';
      } // If MCParticle parent.
    } // if !isRealData
  } // for PFParticle daughters

} // FillPrimaryPFParticle

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void protoana::DUNEPizeroAnaTree::setPiZeroInfo(art::Event const & evt, const std::vector<const simb::MCParticle*>& pi0s) {
  // Clean slate
  // ResetPi0Vars();

  // Loop over all pi0s in the array by index.
  fNMCPi0s = pi0s.size();
  for(unsigned pi0i = 0; pi0i < NMAXPIZEROS && pi0i < fNMCPi0s; ++pi0i) {
    pizero::PiZeroProcess pzproc(*(pi0s[pi0i]), evt, fShowerTag);

    // Set MC object variables
    if(pzproc.allMCSet()) {

      const simb::MCParticle* pi0 = pzproc.pi0();
      const simb::MCParticle* photon1 = pzproc.photon1();
      const simb::MCParticle* photon2 = pzproc.photon2();

      // MC pi0 variables
      fMCPi0ID[pi0i] = pi0->TrackId();
      fMCPi0Energy[pi0i] = pi0->E();
      pi0->Position().Vect().GetXYZ(fMCPi0StartPosition[pi0i]);
      pi0->EndPosition().Vect().GetXYZ(fMCPi0EndPosition[pi0i]);
      fMCPhoton1ID[pi0i] = photon1->TrackId();
      fMCPhoton2ID[pi0i] = photon2->TrackId();
      fMCPhoton1Energy[pi0i] = photon1->E();
      fMCPhoton2Energy[pi0i] = photon2->E();
      photon1->Position().Vect().GetXYZ(fMCPhoton1StartPosition[pi0i]);
      photon2->Position().Vect().GetXYZ(fMCPhoton2StartPosition[pi0i]);
      photon1->EndPosition().Vect().GetXYZ(fMCPhoton1EndPosition[pi0i]);
      photon2->EndPosition().Vect().GetXYZ(fMCPhoton2EndPosition[pi0i]);

    } // if MC set
  } // Loop over pi0s

} // setPiZeroInfo


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void protoana::DUNEPizeroAnaTree::Initialise(){
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  for(int k=0; k < 3; k++){
    fbeamtrackPos[k] = -999.0;
    fbeamtrackEndPos[k] = -999.0;
    fbeamtrackDir[k] = -999.0;
    fprimaryTruth_vtx[k]= -999.0;
    fprimaryVertex[k] = -999.0;
    fprimaryEndPosition[k] = -999.0;
    fprimaryStartPosition[k] = -999.0;
    fprimaryEndDirection[k] = -999.0;
    fprimaryStartDirection[k] = -999.0;
    fprimaryKineticEnergy[k] = -999.0;
    fprimaryRange[k] = -999.0;
  }
  for( int l=0; l<1000; l ++) {
     fbeamtrackPos_at[l][0] = -999.;
     fbeamtrackPos_at[l][1] = -999.;
     fbeamtrackPos_at[l][2] = -999.;
     fbeamtrackMom_at[l][0] = -999.;
     fbeamtrackMom_at[l][1] = -999.;
     fbeamtrackMom_at[l][2] = -999.;
     fbeamtrackMom_at[l][3] = -999.;
  }
  fprimaryShower_nHits =0;
  for( int m=0; m<NMAXHITS; m ++){
     fprimarydEdx[m]= -999.0;
     fprimarydQdx[m]= -999.0;
     fprimary_cal_pos[m][0] = -999.0;
     fprimary_cal_pos[m][1] = -999.0;
     fprimary_cal_pos[m][2] = -999.0;
     fprimary_cal_pitch[m]=-999.0;
     fprimaryResidualRange[m] = -999.0;

     fprimaryShower_hit_w[m] =-999.0;
     fprimaryShower_hit_q[m] =-999.0;
     fprimaryShower_hit_t[m] =-999.0;
     fprimaryShower_hit_pos[m][0] =-999.0;
     fprimaryShower_hit_pos[m][1] =-999.0;
     fprimaryShower_hit_pos[m][2] =-999.0;
  }
  fprimaryTruth_trkID =-999;
  fprimaryTruth_pdg = -999;
  fprimaryTruth_E = -999;
  fprimarynCal = 0;
  fbeamCheckIsMatched = -999;
  fbeamtrigger = -999;
  ftof = -999.0;
  for(int k=0; k < 2; k++){
    fcerenkovPressure[k] = -999.0;
    fcerenkovTime[k] = -999.0;
    fcerenkovStatus[k] = -999;
  }
  fbeamtrackMomentum = -999.0;
  fbeamtrackEnergy = 999.0;
  fbeamtrackPdg = -999;
  fbeamtrackTime = -999.0;
  fbeamtrackID = -999;
  fbeam_ntrjPoints =0;
  fbeamtrackNDaughters= 0;
  fprimaryIstrack = -999;
  fprimaryIsshower = -999;

  fprimaryBDTScore = -999.0;
  fprimaryCNNScore = -999.0;
  fprimaryNHits = 0;
  fprimaryTheta = -999.0;
  fprimaryPhi = -999.0;
  fprimaryLength = -999.0;
  fprimaryMomentum = -999.0;
  fprimaryEndMomentum = -999.0;
  fprimaryOpeningAngle = -999.0;
  fprimaryShowerBestPlane = -999;
  fprimaryShowerEnergy = -999.0;
  fprimaryShowerCharge = -999.0;
  fprimaryShowerMIPEnergy = -999.0;
  fprimaryID = -999;
  fprimaryMomentumByRangeProton = -999.0;
  fprimaryT0 = -999.0;
  fprimaryNDaughters = 0;

  for(unsigned k = 0; k < NMAXDAUGHTERS; ++k) {
    fprimaryDaughterID[k] = -999;
    fprimaryDaughterIstrack[k] = -999;
    fprimaryDaughterIsshower[k] = -999;
    fprimaryDaughterBDTScore[k] = -999.0;
    fprimaryDaughterCNNScore[k] = -999.0;
    fprimaryDaughterNCollHits[k] = 0;
    fprimaryDaughterNHits[k] = 0;
    fprimaryDaughterEnergy[k] = -999.0;
    fprimaryDaughterEnergyFromHits[k] = -999.0;
    fprimaryDaughterMomentum[k] = -999.0;
    fprimaryDaughterEndMomentum[k] = -999.0;
    fprimaryDaughterLength[k] = -999.0;
    for(unsigned l = 0; l < 3; ++l) {
      fprimaryDaughterStartPosition[k][l] = -999.0;
      fprimaryDaughterEndPosition[k][l] = -999.0;
      fprimaryDaughterStartDirection[k][l] = -999.0;
      fprimaryDaughterEndDirection[k][l] = -999.0;
    }
    fprimaryDaughterParentPdg[k] = -999;
    fprimaryDaughterGrandparentPdg[k] = -999;
    fprimaryDaughterParentID[k] = -999;
    fprimaryDaughterGrandparentID[k] = -999;
    fprimaryDaughterParentE[k] = -999.0;
    fprimaryDaughterGrandparentE[k] = -999.0;
    fprimaryDaughterParentStart[k][0] = -999.0;
    fprimaryDaughterParentStart[k][1] = -999.0;
    fprimaryDaughterParentStart[k][2] = -999.0;
    fprimaryDaughterGrandparentStart[k][0] = -999.0;
    fprimaryDaughterGrandparentStart[k][1] = -999.0;
    fprimaryDaughterGrandparentStart[k][2] = -999.0;
    fprimaryDaughterParentEnd[k][0] = -999.0;
    fprimaryDaughterParentEnd[k][1] = -999.0;
    fprimaryDaughterParentEnd[k][2] = -999.0;
    fprimaryDaughterGrandparentEnd[k][0] = -999.0;
    fprimaryDaughterGrandparentEnd[k][1] = -999.0;
    fprimaryDaughterGrandparentEnd[k][2] = -999.0;
    fprimaryDaughterParentEndProcess[k] = -999;
    fprimaryDaughterCompleteness[k] = -999.0;
    fprimaryDaughterPurity[k] = -999.0;
    for(unsigned l = 0; l < NMAXHITS; ++l) {
      // fprimaryDaughterHitdQdx[k][l] = -999.0;
      fprimaryDaughterCaloHitdQdx[k][l] = -999.0;
      fprimaryDaughterCaloHitPitch[k][l] = -999.0;
      fprimaryDaughterOwnHitPitch[k][l] = -999.0;
      fprimaryDaughterHitCharge[k][l] = -999.0;
      // fprimaryDaughterHitdEdxArea[k][l] = -999.0;
      fprimaryDaughterHitWire[k][l] = -999;
      fprimaryDaughterHitTime[k][l] = -999.0;
      fprimaryDaughterHitPosition[k][l][0] = -999.0;
      fprimaryDaughterHitPosition[k][l][1] = -999.0;
      fprimaryDaughterHitPosition[k][l][2] = -999.0;
    }
  }

  // MC pi0 variables.
  for(int pi0i = 0; pi0i < NMAXPIZEROS; ++pi0i) {
    // MC pi0 variables
    fNMCPi0s = -999;
    fMCPi0ID[pi0i] = -999;
    fMCPi0Energy[pi0i] = -999.0;
    for(int i=0; i<3; ++i) {
      fMCPi0StartPosition[pi0i][i] = -999.0;
      fMCPi0EndPosition[pi0i][i] = -999.0;
    }
    fMCPhoton1ID[pi0i] = -999;
    fMCPhoton2ID[pi0i] = -999;
    fMCPhoton1Energy[pi0i] = -999.0;
    fMCPhoton2Energy[pi0i] = -999.0;
    for(int i=0; i<3; ++i) {
      fMCPhoton1StartPosition[pi0i][i] = -999.0;
      fMCPhoton2StartPosition[pi0i][i] = -999.0;
      fMCPhoton1EndPosition[pi0i][i] = -999.0;
      fMCPhoton2EndPosition[pi0i][i] = -999.0;
    }
  }
} // Initialise

DEFINE_ART_MODULE(protoana::DUNEPizeroAnaTree)
