////////////////////////////////////////////////////////////////////////
// Class:  MCS_Tree_Maker
// Module Type: analyzer
// File: MCS_Tree_Maker_module.cc
// Author: Sungbin Oh | sungbino@fnal.gov
// Description: Produe Tree with momentum and segement angles
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"


#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

const int kMaxTracks  = 30;

using namespace std;

namespace dune{

  class MCS_Tree_Maker : public art::EDAnalyzer {
  public:

    explicit MCS_Tree_Maker(fhicl::ParameterSet const& pset);
    virtual ~MCS_Tree_Maker();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    
  private:
    TTree* fEventTree;
    bool     isData;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Double_t evttime; 
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    Int_t    all_trks;
    vector<double>  trackthetaxz;
    vector<double>  trackthetayz;
    vector<double> trkstartx;
    vector<double> trkstarty;
    vector<double> trkstartz;
    vector<double> trkendx;
    vector<double> trkendy;
    vector<double> trkendz;
    vector<double> trklen;
    vector<int>    TrkID;
    vector<vector<int>>    ntrkhits;
    vector<vector<double>> trkdqdx;
    vector<vector<double>>  trkdedx;
    vector<vector<double>>  trkresrange;
    vector<vector<double>>  trkhitx;
    vector<vector<double>>  trkhity;
    vector<vector<double>>  trkhitz;
    vector<vector<double>>  trkpitch;
    
    // == For truth info
    vector<int> true_PID;
    vector<vector<double>> true_momentum;
    vector<vector<double>> true_Px;
    vector<vector<double>> true_Py;
    vector<vector<double>> true_Pz;
    vector<vector<double>> true_hitx;
    vector<vector<double>> true_hity;
    vector<vector<double>> true_hitz;
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;

    protoana::ProtoDUNETruthUtils truthUtil;

  }; 

  //========================================================================
  MCS_Tree_Maker::MCS_Tree_Maker(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false))
  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }
 
  //========================================================================
  MCS_Tree_Maker::~MCS_Tree_Maker(){
  }
  //========================================================================

  //========================================================================
  void MCS_Tree_Maker::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("isData", &isData, "isData/O");
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("year_month_date", &year_month_date,"year_month_date/I");
    fEventTree->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
    fEventTree->Branch("all_trks",&all_trks,"all_trks/I");
    fEventTree->Branch("trackthetaxz","vector<double>", &trackthetaxz);
    fEventTree->Branch("trackthetayz","vector<double>", &trackthetayz);
    fEventTree->Branch("trkstartx","vector<double>", &trkstartx);
    fEventTree->Branch("trkstarty","vector<double>", &trkstarty);
    fEventTree->Branch("trkstartz","vector<double>", &trkstartz);
    fEventTree->Branch("trkendx","vector<double>", &trkendx);
    fEventTree->Branch("trkendy","vector<double>", &trkendy);
    fEventTree->Branch("trkendz","vector<double>", &trkendz);
    fEventTree->Branch("trklen","vector<double>", &trklen);
    fEventTree->Branch("TrkID","vector<int>", &TrkID);
    fEventTree->Branch("ntrkhits","vector<vector<int>>", &ntrkhits);
    fEventTree->Branch("trkdqdx","vector<vector<double>>", &trkdqdx);
    fEventTree->Branch("trkdedx","vector<vector<double>>", &trkdedx);
    fEventTree->Branch("trkresrange","vector<vector<double>>", &trkresrange);
    fEventTree->Branch("trkhitx","vector<vector<double>>", &trkhitx);
    fEventTree->Branch("trkhity","vector<vector<double>>", &trkhity);
    fEventTree->Branch("trkhitz","vector<vector<double>>", &trkhitz);
    fEventTree->Branch("trkpitch","vector<vector<double>>", &trkpitch);

    // == For truth info
    fEventTree->Branch("true_PID","vector<int>",&true_PID);
    fEventTree->Branch("true_momentum","vector<vector<double>>", &true_momentum); 
    fEventTree->Branch("true_Px","vector<vector<double>>", &true_Px);
    fEventTree->Branch("true_Py","vector<vector<double>>", &true_Py);
    fEventTree->Branch("true_Pz","vector<vector<double>>", &true_Pz);
    fEventTree->Branch("true_hitx","vector<vector<double>>", &true_hitx);
    fEventTree->Branch("true_hity","vector<vector<double>>", &true_hity);
    fEventTree->Branch("true_hitz","vector<vector<double>>", &true_hitz);
  }

  //========================================================================
  void MCS_Tree_Maker::endJob(){     

  }

  //========================================================================
  void MCS_Tree_Maker::beginRun(const art::Run&){
    mf::LogInfo("MCS_Tree_Maker")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void MCS_Tree_Maker::analyze( const art::Event& evt){
    reset();  
    std::vector<art::Ptr<recob::Track> > tracklist;
    auto trackListHandle = evt.getHandle< std::vector<recob::Track> >("pandoraTrack");
    if (trackListHandle){
      art::fill_ptr_vector(tracklist, trackListHandle);
    }
    else return;

    std::vector<art::Ptr<recob::PFParticle> > pfplist;    
    auto PFPListHandle = evt.getHandle< std::vector<recob::PFParticle> >("pandora"); 
    if (PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);
  
    std::vector<art::Ptr<recob::Hit>> hitlist;
    auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel) ; // to get information about the hits
    if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);
       
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits

    art::FindManyP<recob::Hit> fmtht(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
    art::FindManyP<recob::Track> thass(hitListHandle, evt, fTrackModuleLabel); //to associate hit just trying
    
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    if (!fmcal.isValid()){
      std::cout<<"art::Assns<recob::Track,anab::Calorimetry,void> with module label: "<<fCalorimetryModuleLabel<<" not found."<<std::endl;
      return;
    }

    art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, evt ,"pandora");
   
  
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,evt,"pandoraTrack");
    art::FindManyP<anab::T0> fmT0(trackListHandle, evt ,"pmtrack");

    // Get services
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    isData = evt.isRealData();
    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime=tts.AsDouble();
     
    UInt_t year=0;
    UInt_t month=0;
    UInt_t day=0;
     
    year_month_date=tts.GetDate(kTRUE,0,&year,&month,&day);
     
    UInt_t hour=0;
    UInt_t min=0;
    UInt_t sec=0;
     
    hour_min_sec=tts.GetTime(kTRUE,0,&hour,&min,&sec);
  
    all_trks=0;
    // size_t NHits=hitlist.size();
    size_t NTracks = tracklist.size();
    for(size_t i=0; i<NTracks;++i){
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      if(fTrackModuleLabel=="pandoraTrack"){
	std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(i);
	if(!pfps.size()) continue;
	std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
	if(!t0s.size()) continue;
	//auto t0 = t0s.at(0);
	// double t_zero=t0->Time();
      }
       
      if(fTrackModuleLabel=="pmtrack"){
	std::vector<art::Ptr<anab::T0>> T0s=fmT0.at(i);
	if(T0s.size()==0)
	  continue;
      }
      all_trks++;

      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i); // == calos : dq/dx, dE/dx, hit positions and pitches
      const recob::Track& track = *ptrack;
      auto pos = track.Vertex();
      auto dir_start = track.VertexDirection();
      auto dir_end   = track.EndDirection();
      auto end = track.End();
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());  
      float startx=pos.X();   
      float starty=pos.Y();
      float startz=pos.Z();
      float endx=end.X();
      float endy=end.Y();
      float endz=end.Z();
      float tracklength=track.Length();
      size_t count = 0;

      /*************************************filling the values****************************/
      trackthetaxz.push_back(theta_xz);
      trackthetayz.push_back(theta_yz);
      trkstartx.push_back(startx);
      trkstarty.push_back(starty);
      trkstartz.push_back(startz);
      trkendx.push_back(endx);
      trkendy.push_back(endy);
      trkendz.push_back(endz);
      trklen.push_back(tracklength);
      TrkID.push_back(track.ID());
      for(size_t ical = 0; ical<calos.size(); ++ical){
	if(!calos[ical]) continue;
	if(!calos[ical]->PlaneID().isValid) continue;
	int planenum = calos[ical]->PlaneID().Plane;
	if(planenum != 2) continue; // == Save only collection view hits
	const size_t NHits = calos[ical] -> dEdx().size();
	ntrkhits[cross_trks-1][planenum]=int(NHits);
	vector<double> this_trkdqdx;
	vector<double> this_trkdedx;
	vector<double> this_trkresrange;
	vector<double> this_trkhitx;
	vector<double> this_trkhity;
	vector<double> this_trkhitz;
	vector<double> this_trkpitch;
	for(size_t iHit = 0; iHit < NHits; ++iHit){
	  const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
	  this_trkdqdx.push_back((calos[ical] -> dQdx())[iHit]);
	  this_trkdedx.push_back((calos[ical] -> dEdx())[iHit]);
	  this_trkresrange.push_back((calos[ical]->ResidualRange())[iHit]);
	  this_trkhitx.push_back(TrkPos.X());
	  this_trkhity.push_back(TrkPos.Y());
          this_trkhitz.push_back(TrkPos.Z());
	  this_trkpitch.push_back((calos[ical]->TrkPitchVec())[iHit]);
	} // loop over iHit..

	trkdqdx.push_back(this_trkdqdx);
	trkdedx.push_back(this_trkdedx);
	trkresrange.push_back(this_trkresrange);
	trkhitx.push_back(this_trkhitx);
	trkhity.push_back(this_trkhity);
	trkhitz.push_back(this_trkhitz);
	trkpitch.push_back(this_trkpitch);

      } // loop over ical 2nd time...
      if (!isData){

        auto particles = truthUtil.GetMCParticleListFromRecoTrack(clockData, *tracklist[i], evt, fTrackModuleLabel);

        for (auto & part : particles){
          //cout<<(part.first)->PdgCode()<<" "<<part.second<<endl;
          if (!part.first) continue;
	  true_PID.push_back(((part.first)->PdgCode()));
	
	  vector<double> this_true_momentum;
	  vector<double> this_true_Px;
	  vector<double> this_true_Py;
          vector<double> this_true_Pz;
	  vector<double> this_true_hitx;
	  vector<double> this_true_hity;
	  vector<double> this_true_hitz;
  
	  int N_positions = (part.first)->Position().size();
	  const simb::MCTrajectory & true_trajectory = (part.first)->Trajectory();
	  auto true_proc_map = true_trajectory.TrajectoryProcesses();
	  for( auto itProc = true_proc_map.begin(); itProc != true_proc_map.end(); ++itProc ){
	    int index = itProc->first;
	    //std::string process = true_trajectory.KeyToProcess(itProc->second);

	    double this_hit_Px = true_beam_trajectory.Px(index);
	    double this_hit_Py = true_beam_trajectory.Py(index);
            double this_hit_Pz = true_beam_trajectory.Pz(index);
	    double this_hit_P = sqrt(this_hit_Px*this_hit_Px + this_hit_Py*this_hit_Py + this_hit_Pz*this_hit_Pz);
	    this_true_momentum.push_back(this_hit_P);
	    this_true_Px.push_back(this_hit_Px);
	    this_true_Py.push_back(this_hit_Py);
            this_true_Pz.push_back(this_hit_Pz);
	 
	    this_true_hitx.push_back(true_beam_trajectory.X(index));
            this_true_hity.push_back(true_beam_trajectory.Y(index));
            this_true_hitz.push_back(true_beam_trajectory.Z(index));
	  }

	  true_momentum.push_back(this_true_momentum);
	  true_Px.push_back(this_true_Px);
          true_Py.push_back(this_true_Py);
          true_Pz.push_back(this_true_Pz);
	  true_hitx.push_back(this_true_hitx);
	  true_hity.push_back(this_true_hity);
          true_hitz.push_back(this_true_hitz);

          break;
        }
      }
    } // loop over trks...
  
    fEventTree->Fill();
  } // end of analyze function
     
    /////////////////// Defintion of reset function ///////////
  void MCS_Tree_Maker::reset(){
    isData = false;
    run = -99999;
    subrun = -99999;
    event = -99999;
    evttime = -99999;
    cross_trks = -99999;
    all_trks = -99999;
    year_month_date=-99999;
    hour_min_sec=-99999;

    trackthetaxz.clear();
    trackthetayz.clear();
    trkstartx.clear();
    trkstarty.clear();
    trkstartz.clear();
    trkendx.clear();
    trkendy.clear();
    trkendz.clear();
    trklen.clear();
    TrkID.clear();
    ntrkhits.clear();
    trkdqdx.clear();
    trkdedx.clear();
    trkresrange.clear();
    trkhitx.clear();
    trkhity.clear();
    trkhitz.clear();
    trkpitch.clear();

    // == For truth info
    true_PID.clear();
    true_momentum.clear();
    true_Px.clear();
    true_Py.clear();
    true_Pz.clear();
    true_hitx.clear();
    true_hity.clear();
    true_hitz.clear();
  }
  //////////////////////// End of definition ///////////////
    
  DEFINE_ART_MODULE(MCS_Tree_Maker)
}
