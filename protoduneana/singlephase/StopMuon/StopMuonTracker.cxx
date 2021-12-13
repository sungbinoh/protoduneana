////////////////////////////////////////////////////////////////////////
// Class:       StopMuonTracker
// Plugin Type: analyzer (art v3_02_06)
// File:        StopMuonTracker.cxx
//
// Generated at Thu Feb 27 20:32:44 2020 by Jason Stock using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "StopMuonTracker.h"

StopMuonTracker::StopMuonTracker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPFPModuleLabel(p.get<art::InputTag>("PFPModuleLabel", "default")),
  fTrackModuleLabel(p.get<art::InputTag>("TrackModuleLabel", "default")),
  fOpHitModuleLabel(p.get<art::InputTag>("OpHitModuleLabel", "default")),
  fFlashModuleLabel(p.get<art::InputTag>("FlashModuleLabel", "default")),
  fOpHitTriggerLabel(p.get<art::InputTag>("OpHitTriggerModuleLabel", "default")),
  fFlashTriggerLabel(p.get<art::InputTag>("FlashTriggerModuleLabel", "default")),
  fSpacePointSolverModuleLabel(p.get<art::InputTag>("SpacePointLabel", "default")),
  fHitModuleLabel(p.get<art::InputTag>("HitModuleLabel", "default")),
  fShowerModuleLabel(p.get<art::InputTag>("ShowerModuleLabel", "default"))

  // More initializers here.
{
  std::cout<<"Printing Module Labels:\n"
    <<"PFP: "<<fPFPModuleLabel.label()<<"\n"
    <<"Track: "<<fTrackModuleLabel.label()<<"\n"
    <<"OpHit: "<<fOpHitModuleLabel.label()<<"\n"
    <<"Flash: "<<fFlashModuleLabel.label()<<"\n"
    <<"OpHitTrigger: "<<fOpHitTriggerLabel.label()<<"\n"
    <<"FlashTrigger: "<<fFlashTriggerLabel.label()<<"\n"
    <<"SpacePoint: "<<fSpacePointSolverModuleLabel.label()<<"\n"
    <<"Hit: "<<fHitModuleLabel.label()<<"\n"
    <<"Shower: "<<fShowerModuleLabel.label()<<"\n";
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

  StopMuonTracker::StopMuonTracker( art::EDAnalyzer::Table<StopMuonTrackerConfig> const& config)
: art::EDAnalyzer(config),
  fPFPModuleLabel(config().PFPModuleLabel()),
  fTrackModuleLabel(config().TrackModuleLabel()),
  fOpHitModuleLabel(config().OpHitModuleLabel()),
  fFlashModuleLabel(config().FlashModuleLabel()),
  fOpHitTriggerLabel(config().OpHitTriggerLabel()),
  fFlashTriggerLabel(config().FlashTriggerLabel()),
  fSpacePointSolverModuleLabel(config().SpacePointLabel()),
  fHitModuleLabel(config().HitModuleLabel()),
  fShowerModuleLabel(config().ShowerModuleLabel())
{
  std::cout<<"Printing Module Labels:\n"
    <<"PFP: "<<fPFPModuleLabel.label()<<"\n"
    <<"Track: "<<fTrackModuleLabel.label()<<"\n"
    <<"OpHit: "<<fOpHitModuleLabel.label()<<"\n"
    <<"Flash: "<<fFlashModuleLabel.label()<<"\n"
    <<"OpHitTrigger: "<<fOpHitTriggerLabel.label()<<"\n"
    <<"FlashTrigger: "<<fFlashTriggerLabel.label()<<"\n"
    <<"SpacePoint: "<<fSpacePointSolverModuleLabel.label()<<"\n"
    <<"Hit: "<<fHitModuleLabel.label()<<"\n"
    <<"Shower: "<<fShowerModuleLabel.label()<<"\n";
}

void StopMuonTracker::analyze(art::Event const& e)
{

  std::cout<<"EVENT_NOTE: Starting"<<"\n";
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);


  Int_t skipped_events_nf         =-10;
  Int_t skipped_events_oh         =-10;
  Int_t skipped_events_nt         =-10;
  Int_t track_has_pfp             =-10;
  Int_t track_has_t0              =-10;
  Int_t track_in_bounds           =-10;
  Int_t track_crosses_cathode     =-10;
  Int_t track_passes_fiducial_cut =-10;
  Int_t track_not_broken          =-10;
  Int_t track_not_short           =-10;
  Int_t track_has_hits            =-10;
  Int_t track_has_no_pfp          =-10;
  Int_t track_has_no_t0           =-10;
  Int_t track_not_in_bounds       =-10;
  Int_t track_no_crosses_cathode  =-10;
  Int_t track_fails_fiducial_cut  =-10;
  Int_t track_is_broken           =-10;
  Int_t track_short               =-10;
  Int_t track_has_no_hits         =-10;
  Int_t track_has_no_ohits        =-10;
  Int_t track_has_no_oflash       =-10;
  Int_t track_has_ohits           =-10;
  Int_t track_has_opflash         =-10;


  priv_TrackBranchBuffer.clear();

  priv_TrackBranchBuffer.set_run(e.id().run());
  priv_TrackBranchBuffer.set_subrun(e.id().subRun());
  priv_TrackBranchBuffer.set_event(e.id().event());
  {
    TTimeStamp tmp(e.time().timeHigh(), e.time().timeLow());
    priv_TrackBranchBuffer.set_time(tmp.AsDouble());
  }

  std::cout<<"EVENT_NOTE: Logged run,subrun,event\n";
  art::Handle< std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  art::Handle< std::vector<recob::Track >> trackListHandle;
  std::vector<art::Ptr<recob::Track>> tracklist;
  art::Handle< std::vector<recob::PFParticle>> PFPListHandle;
  std::vector<art::Ptr<recob::PFParticle>>pfplist;
  art::Handle< std::vector<recob::OpHit> > ophitListHandle;
  std::vector<art::Ptr<recob::OpHit> > ophitList;
  art::Handle< std::vector<recob::OpFlash> > flashListHandle;
  std::vector<art::Ptr<recob::OpFlash> > flashList;
  art::Handle< std::vector<recob::OpFlash>> flashTriggerHandle;
  art::Handle< std::vector<recob::OpHit>>   ophitTriggerHandle;
  if(e.getByLabel(fPFPModuleLabel,PFPListHandle) && e.getByLabel(fTrackModuleLabel, trackListHandle) && e.getByLabel(fOpHitModuleLabel, ophitListHandle) && e.getByLabel(fFlashModuleLabel, flashListHandle) && e.getByLabel(fHitModuleLabel, hitListHandle) && e.getByLabel(fShowerModuleLabel, showerListHandle ))
  {
    art::fill_ptr_vector(tracklist, trackListHandle);
    art::fill_ptr_vector(pfplist, PFPListHandle);
    art::fill_ptr_vector(ophitList, ophitListHandle);
    art::fill_ptr_vector(flashList, flashListHandle);
    art::fill_ptr_vector(hitlist, hitListHandle);
    art::fill_ptr_vector(showerlist, showerListHandle);
  }
  std::cout<<"EVENT_NOTE: Handle set\n";
  double flashTrigT(-9999999.9), ophitTrigT(-9999999.9);

  if( !e.getByLabel(fFlashTriggerLabel, flashTriggerHandle) ){
    std::cout<<"FlashTriggerLabel: "<<fFlashTriggerLabel.label()<<"\n";
    std::cout<<"Failed to get Flash Trigger. Skipping Event!\n";
    //std::string tmpst=(std::to_string(e.id().run())+std::to_string(e.id().subRun())+std::to_string(e.id().event()));
    //skipped_events_nf=std::stoi(tmpst);
    skipped_events_nf=1;
    priv_skipped_events_nf_Buffer         = skipped_events_nf        ;
    priv_skipped_events_oh_Buffer         = skipped_events_oh        ;
    priv_skipped_events_nt_Buffer         = skipped_events_nt        ;
    priv_track_has_pfp_Buffer             = track_has_pfp            ;
    priv_track_has_t0_Buffer              = track_has_t0             ;
    priv_track_in_bounds_Buffer           = track_in_bounds          ;
    priv_track_crosses_cathode_Buffer     = track_crosses_cathode    ;
    priv_track_passes_fiducial_cut_Buffer = track_passes_fiducial_cut;
    priv_track_not_broken_Buffer          = track_not_broken         ;
    priv_track_not_short_Buffer           = track_not_short          ;
    priv_track_has_hits_Buffer            = track_has_hits           ;
    priv_track_has_opflash_Buffer         = track_has_opflash        ;
    priv_track_has_no_pfp_Buffer          = track_has_no_pfp         ;
    priv_track_has_no_t0_Buffer           = track_has_no_t0          ;
    priv_track_not_in_bounds_Buffer       = track_not_in_bounds      ;
    priv_track_no_crosses_cathode_Buffer  = track_no_crosses_cathode ;
    priv_track_fails_fiducial_cut_Buffer  = track_fails_fiducial_cut ;
    priv_track_is_broken_Buffer           = track_is_broken          ;
    priv_track_short_Buffer               = track_short              ;
    priv_track_has_no_hits_Buffer         = track_has_no_hits        ;
    priv_track_has_no_ohits_Buffer        = track_has_no_ohits       ;
    priv_track_has_no_oflash_Buffer       = track_has_no_oflash      ;
    priv_track_has_ohits_Buffer           = track_has_ohits          ;
    priv_StatTree->Fill(); 
    return;
  }
  if( !e.getByLabel(fOpHitTriggerLabel, ophitTriggerHandle))
  {
    //std::string tmpst=(std::to_string(e.id().run())+std::to_string(e.id().subRun())+std::to_string(e.id().event()));
    //skipped_events_oh=std::stoi(tmpst);
    skipped_events_oh=1;
    std::cout<<"Skipping event! No ophits T0\n";
    priv_skipped_events_nf_Buffer         = skipped_events_nf        ;
    priv_skipped_events_oh_Buffer         = skipped_events_oh        ;
    priv_skipped_events_nt_Buffer         = skipped_events_nt        ;
    priv_track_has_pfp_Buffer             = track_has_pfp            ;
    priv_track_has_t0_Buffer              = track_has_t0             ;
    priv_track_in_bounds_Buffer           = track_in_bounds          ;
    priv_track_crosses_cathode_Buffer     = track_crosses_cathode    ;
    priv_track_passes_fiducial_cut_Buffer = track_passes_fiducial_cut;
    priv_track_not_broken_Buffer          = track_not_broken         ;
    priv_track_not_short_Buffer           = track_not_short          ;
    priv_track_has_hits_Buffer            = track_has_hits           ;
    priv_track_has_opflash_Buffer         = track_has_opflash        ;
    priv_track_has_no_pfp_Buffer          = track_has_no_pfp         ;
    priv_track_has_no_t0_Buffer           = track_has_no_t0          ;
    priv_track_not_in_bounds_Buffer       = track_not_in_bounds      ;
    priv_track_no_crosses_cathode_Buffer  = track_no_crosses_cathode ;
    priv_track_fails_fiducial_cut_Buffer  = track_fails_fiducial_cut ;
    priv_track_is_broken_Buffer           = track_is_broken          ;
    priv_track_short_Buffer               = track_short              ;
    priv_track_has_no_hits_Buffer         = track_has_no_hits        ;
    priv_track_has_no_ohits_Buffer        = track_has_no_ohits       ;
    priv_track_has_no_oflash_Buffer       = track_has_no_oflash      ;
    priv_track_has_ohits_Buffer           = track_has_ohits          ;
    priv_StatTree->Fill(); 
    return;
  } 
  if( flashTriggerHandle->empty() || ophitTriggerHandle->empty() )
  {
    std::cout<<"Skipping Event. Flash or Trigger is missing.\n";
    return;
  }
  flashTrigT = flashTriggerHandle->at(0).Time() ;
  ophitTrigT = ophitTriggerHandle->at(0).PeakTime() ;
  if(dbg)
  {
    std::cout<<"FlashT0 - OpHitT0 = "<<flashTrigT-ophitTrigT<<"\n";
    std::cout<<"FlashT0: "<<flashTrigT<<" OpHitT0: "<<ophitTrigT<<"\n";
  }
  if(flashTrigT < -9999999.0 || ophitTrigT < -9999999.0)
  {
    //std::string tmpst=(std::to_string(e.id().run())+std::to_string(e.id().subRun())+std::to_string(e.id().event()));
    //skipped_events_nt=std::stoi(tmpst);
    skipped_events_nt=1;
    std::cout<<"Skipping event! Missing optical T0\n";
    priv_skipped_events_nf_Buffer         = skipped_events_nf        ;
    priv_skipped_events_oh_Buffer         = skipped_events_oh        ;
    priv_skipped_events_nt_Buffer         = skipped_events_nt        ;
    priv_track_has_pfp_Buffer             = track_has_pfp            ;
    priv_track_has_t0_Buffer              = track_has_t0             ;
    priv_track_in_bounds_Buffer           = track_in_bounds          ;
    priv_track_crosses_cathode_Buffer     = track_crosses_cathode    ;
    priv_track_passes_fiducial_cut_Buffer = track_passes_fiducial_cut;
    priv_track_not_broken_Buffer          = track_not_broken         ;
    priv_track_not_short_Buffer           = track_not_short          ;
    priv_track_has_hits_Buffer            = track_has_hits           ;
    priv_track_has_opflash_Buffer         = track_has_opflash        ;
    priv_track_has_no_pfp_Buffer          = track_has_no_pfp         ;
    priv_track_has_no_t0_Buffer           = track_has_no_t0          ;
    priv_track_not_in_bounds_Buffer       = track_not_in_bounds      ;
    priv_track_no_crosses_cathode_Buffer  = track_no_crosses_cathode ;
    priv_track_fails_fiducial_cut_Buffer  = track_fails_fiducial_cut ;
    priv_track_is_broken_Buffer           = track_is_broken          ;
    priv_track_short_Buffer               = track_short              ;
    priv_track_has_no_hits_Buffer         = track_has_no_hits        ;
    priv_track_has_no_ohits_Buffer        = track_has_no_ohits       ;
    priv_track_has_no_oflash_Buffer       = track_has_no_oflash      ;
    priv_track_has_ohits_Buffer           = track_has_ohits          ;
    priv_StatTree->Fill(); 
    return;
    // Put a proper mf log here. Warning level;
  }


  art::FindManyP<recob::Hit> fmp_hit (trackListHandle, e, fTrackModuleLabel);
  art::FindManyP<recob::PFParticle> fmp_pfp (trackListHandle, e, fTrackModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_thm(trackListHandle, e, fTrackModuleLabel);
  art::FindManyP<recob::SpacePoint> fmpsp(hitListHandle, e, fSpacePointSolverModuleLabel);
  art::FindManyP<recob::PFParticle> pfp_shwr_assn (showerListHandle, e, fShowerModuleLabel);
  art::FindManyP<anab::T0> shwr_t0_assn_v (PFPListHandle, e, fPFPModuleLabel);

  std::cout<<"Starting Track Loop.\n";
  for(size_t i=0; i<tracklist.size(); i++)
  {
    /* art::Ptr<recob::Track> ptrack(trackListHandle, i); */
    art::Ptr<recob::Track> ptrack = tracklist.at(i);
    const recob::Track& track = *ptrack;
    std::vector<art::Ptr<recob::PFParticle>> pfps=fmp_pfp.at(i);

    track_has_no_pfp=(track_has_no_pfp<0)?0:track_has_no_pfp;
    track_has_pfp=(track_has_pfp<0)?0:track_has_pfp;
    if(pfps.empty()){
      ++track_has_no_pfp;
      if(dbg)
      continue; //No flow, no T0. Skip
    }
    ++track_has_pfp;
    art::FindManyP<anab::T0> fmp_t0(PFPListHandle, e, fPFPModuleLabel);
    std::vector<art::Ptr<anab::T0>> t0s = fmp_t0.at(pfps[0].key());
    track_has_no_t0=(track_has_no_t0<0)?0:track_has_no_t0;
    track_has_t0=(track_has_t0<0)?0:track_has_t0;
    if(t0s.empty()){
      ++track_has_no_t0;
      continue; //No T0. Skip
    }
    ++track_has_t0;
    //Include units here.
    double T00=double(t0s.at(0)->Time());
    if(dbg){std::cout<<"T00: "<<T00<<"\n";}
    TVector3 dr_s = track.DirectionAtPoint<TVector3>(track.FirstValidPoint());
    TVector3 dr_e = track.DirectionAtPoint<TVector3>(track.LastValidPoint());
    TVector3 ps_s(track.LocationAtPoint(track.FirstValidPoint()).X(), track.LocationAtPoint(track.FirstValidPoint()).Y(), track.LocationAtPoint(track.FirstValidPoint()).Z()); //inherited from Aleena
    TVector3 ps_e(track.LocationAtPoint(track.LastValidPoint()).X(), track.LocationAtPoint(track.LastValidPoint()).Y(), track.LocationAtPoint(track.LastValidPoint()).Z()); //inherited from Aleena
//    TVector3 ps_e = track.End<TVector3>();

    track_not_in_bounds=(track_not_in_bounds<0)?0:track_not_in_bounds;
    track_in_bounds=(track_in_bounds<0)?0:track_in_bounds;
    if( !MFL::IsPointInBounds(ps_s) && !MFL::IsPointInBounds(ps_e) ) {
      ++track_not_in_bounds;
      continue;
    }
    ++track_in_bounds;
    MFL::ForceDirection(ps_s,ps_e,dr_s,dr_e); //Based on SwitchEndPoints from Aleena
    //First, ps_s.X*ps_e.X<0, second, endZ outside APA{neg,pos}bound{1,2}
    track_no_crosses_cathode=(track_no_crosses_cathode<0)?0:track_no_crosses_cathode;
    track_crosses_cathode=(track_crosses_cathode<0)?0:track_crosses_cathode;
    if( !MFL::CathCross(ps_s,ps_e) ){
      ++track_no_crosses_cathode;
      continue;
    }
    ++track_crosses_cathode;
    track_fails_fiducial_cut=(track_fails_fiducial_cut<0)?0:track_fails_fiducial_cut;
    track_passes_fiducial_cut=(track_passes_fiducial_cut<0)?0:track_passes_fiducial_cut;
    if( !MFL::FidCut(ps_e) ){
      ++track_fails_fiducial_cut;
      continue; //FV, EndY, and EndZ cuts all included.
    }
    ++track_passes_fiducial_cut;

    track_is_broken=(track_is_broken<0)?0:track_is_broken;
    track_not_broken=(track_not_broken<0)?0:track_not_broken;
    if( MFL::BrokenTrack(ptrack, tracklist) ) {
      ++track_is_broken;
      continue;
    }
    ++track_not_broken;

    track_short=(track_short<0)?0:track_short;
    track_not_short=(track_not_short<0)?0:track_not_short;
    if( track.Length() < trk_len_min_cut ){
      ++track_short;
      continue; //Track too short.
    }
    ++track_not_short;

    track_has_no_hits=(track_has_no_hits<0)?0:track_has_no_hits;
    track_has_hits=(track_has_hits<0)?0:track_has_hits;
    std::vector<ana::hit> hits = MFL::GetHitsForTrack(ptrack, e, fTrackModuleLabel); //HitTimesCut is in GetHitsForTrack. Hits sort in GetHitsForTrack.
    if(hits.empty()){ 
      ++track_has_no_hits;
      continue; 
    }
    ++track_has_hits;
    

    track_has_no_ohits=(track_has_no_ohits<0)?0:track_has_no_ohits;
    track_has_ohits=(track_has_ohits<0)?0:track_has_ohits;
    std::vector<ana::ohit> ohits = MFL::GetOpHitsForTrack(ptrack, T00, ophitTrigT, ophitList);
    if(ohits.empty() ){
      ++track_has_no_ohits;
      continue;
    }
    ++track_has_ohits;

    track_has_no_oflash=(track_has_no_oflash<0)?0:track_has_no_oflash;
    track_has_opflash=(track_has_opflash<0)?0:track_has_opflash;
    std::vector<ana::oflash> oflashs = MFL::GetOpFlashesForTrack(ptrack, T00, flashTrigT, flashList);
    if( oflashs.empty() ){
      ++track_has_no_oflash;
      continue;
    }
    ++track_has_opflash;

    ana::smtrec tmpr(i, std::make_pair(MFL::MinTime(hits), T00/1000.0), false, ana::smtrec::track_type_t::stop_muon, track.StartMomentum()  );
    tmpr.add_trigger(T00);
    for(size_t d=0; d<track.NPoints(); d++)
    { 
      if(track.HasValidPoint(d))
      {
        auto & sp = track.LocationAtPoint(d);
        auto mom = sqrt(track.MomentumVectorAtPoint(d).Mag2());
        tmpr.add_point(sp.X(), sp.Y(), sp.Z(), mom); 
      }
    }
    for(auto& hit: hits)
    { tmpr.add_hit(hit);}
    for(auto& ophit : ohits)
    { tmpr.add_ophit(ophit); }
    for(auto &flash : oflashs )
    { tmpr.add_opflash(flash); }
    priv_TrackBranchBuffer.add_smt(tmpr);


    std::cout<<"Filling Trees.\n";
    priv_TrackTree->Fill();
    //Write out Muon info here
    //What do I want to write out? 
    //Muon. Space points (all).
    //Start.
    //Stop.
    //Flashes.
    //Ophits.
    //Hits.
    //Write it all out in a common clock for each event.

    //Tag muons with Michels.
    //HERE TODO: Do I need a fmp_thm isvalid check?
    //I don't want to do this for now. I can survive without michel tagging.
    //std::vector<std::vector<double>> hasMichel = MFL::HasMichel( tracklist, fmp_thm, ps_s, ps_e, i, T00, fmpsp, trackListHandle, hitlist, showerlist, pfp_shwr_assn, shwr_t0_assn_v, );

    //Check for michels, and tag weather track has one or not for write

    //Store this track.
  }
  priv_track_has_pfp_Buffer             = track_has_pfp            ;
  priv_track_has_t0_Buffer              = track_has_t0             ;
  priv_track_in_bounds_Buffer           = track_in_bounds          ;
  priv_track_crosses_cathode_Buffer     = track_crosses_cathode    ;
  priv_track_passes_fiducial_cut_Buffer = track_passes_fiducial_cut;
  priv_track_not_broken_Buffer          = track_not_broken         ;
  priv_track_not_short_Buffer           = track_not_short          ;
  priv_track_has_hits_Buffer            = track_has_hits           ;
  priv_track_has_opflash_Buffer         = track_has_opflash        ;
  priv_track_has_no_pfp_Buffer          = track_has_no_pfp         ;
  priv_track_has_no_t0_Buffer           = track_has_no_t0          ;
  priv_track_not_in_bounds_Buffer       = track_not_in_bounds      ;
  priv_track_no_crosses_cathode_Buffer  = track_no_crosses_cathode ;
  priv_track_fails_fiducial_cut_Buffer  = track_fails_fiducial_cut ;
  priv_track_is_broken_Buffer           = track_is_broken          ;
  priv_track_short_Buffer               = track_short              ;
  priv_track_has_no_hits_Buffer         = track_has_no_hits        ;
  priv_track_has_no_ohits_Buffer        = track_has_no_ohits       ;
  priv_track_has_no_oflash_Buffer       = track_has_no_oflash      ;
  priv_track_has_ohits_Buffer           = track_has_ohits          ;
  priv_StatTree->Fill(); 
  //write out
}

void StopMuonTracker::beginJob()
{
  art::ServiceHandle<art::TFileService> priv_tfs;
  priv_TrackTree = priv_tfs->make<TTree>("TrackRecords",
      "A tree recording muon tracks");
  priv_StatTree = priv_tfs->make<TTree>("StatRecords",
      "A tree recording events stats");
  priv_TrackBranch = priv_TrackTree->Branch("TrackBranch",
      &priv_TrackBranchBuffer, 8000, 99);
  priv_StatBranch1  = priv_StatTree->Branch("SkippedEventsNoFlash",
      &priv_skipped_events_nf_Buffer );
  priv_StatBranch1  = priv_StatTree->Branch("SkippedEventsNoOpHits",
      &priv_skipped_events_oh_Buffer );
  priv_StatBranch2  = priv_StatTree->Branch("SkippedEventsNoT0",
      &priv_skipped_events_nt_Buffer );
  priv_StatBranch3  = priv_StatTree->Branch("TrackHasPFP",
      &priv_track_has_pfp_Buffer );
  priv_StatBranch4  = priv_StatTree->Branch("TrackHasT0",
      &priv_track_has_t0_Buffer );
  priv_StatBranch5  = priv_StatTree->Branch("TrackInBounds",
      &priv_track_in_bounds_Buffer );
  priv_StatBranch6  = priv_StatTree->Branch("TrackCrossesCathode",
      &priv_track_crosses_cathode_Buffer );
  priv_StatBranch7  = priv_StatTree->Branch("TrackPassesFiducialCut",
      &priv_track_passes_fiducial_cut_Buffer );
  priv_StatBranch8  = priv_StatTree->Branch("TrackNotBroken",
      &priv_track_not_broken_Buffer );
  priv_StatBranch9  = priv_StatTree->Branch("TrackNotShort",
      &priv_track_not_short_Buffer );
  priv_StatBranch10 = priv_StatTree->Branch("TrackHasHits",
      &priv_track_has_hits_Buffer );
  priv_StatBranch10 = priv_StatTree->Branch("TrackHasOpFlash",
      &priv_track_has_opflash_Buffer );
  priv_StatBranch11 = priv_StatTree->Branch("TrackHasNoPFP",
      &priv_track_has_no_pfp_Buffer);
  priv_StatBranch12 = priv_StatTree->Branch("TrackHasNoT0",
      &priv_track_has_no_t0_Buffer);
  priv_StatBranch13 = priv_StatTree->Branch("TrackOutOfBounds",
      &priv_track_not_in_bounds_Buffer);
  priv_StatBranch14 = priv_StatTree->Branch("TrackNoCrossCathode",
      &priv_track_no_crosses_cathode_Buffer);
  priv_StatBranch15 = priv_StatTree->Branch("TrackFailsFiducialCut",
      &priv_track_fails_fiducial_cut_Buffer);
  priv_StatBranch16 = priv_StatTree->Branch("TrackIsBroken",
      &priv_track_is_broken_Buffer);
  priv_StatBranch17 = priv_StatTree->Branch("TrackTooShort",
      &priv_track_short_Buffer);
  priv_StatBranch18 = priv_StatTree->Branch("TrackHasNoHits",
      &priv_track_has_no_hits_Buffer);
  priv_StatBranch19 = priv_StatTree->Branch("TrackHasNoOpHits",
      &priv_track_has_no_ohits_Buffer);
  priv_StatBranch20 = priv_StatTree->Branch("TrackHasNoFlash",
      &priv_track_has_no_oflash_Buffer);
  priv_StatBranch21 = priv_StatTree->Branch("TrackHasOpHits",
      &priv_track_has_ohits_Buffer);
}

void StopMuonTracker::beginRun(art::Run const& r)
{
  // Implementation of optional member function here.
}

void StopMuonTracker::endJob()
{
  // Implementation of optional member function here.
}

void StopMuonTracker::endRun(art::Run const& r)
{
  // Implementation of optional member function here.
}


