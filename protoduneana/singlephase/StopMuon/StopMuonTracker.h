#ifndef ___StopMuonTracker_h
#define ___StopMuonTracker_h
////////////////////////////////////////////////////////////////////////
// Class:       StopMuonTracker
// Plugin Type: analyzer (art v3_02_06)
// File:        StopMuonTracker.h
//
// Generated at Thu Feb 27 20:32:44 2020 by Jason Stock using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "dune/DuneObj/SMTracks.h"

#include <TTree.h>
#include <TBranch.h>
#include <TTimeStamp.h>

#include "MuonFinderLib.h"

class StopMuonTracker : public art::EDAnalyzer {
  public:

    struct StopMuonTrackerConfig{
      fhicl::Atom<art::InputTag> PFPModuleLabel{
        fhicl::Name("PFPModuleLabel"),
          fhicl::Comment("Module that created PFPs"),
          "default"};
      fhicl::Atom<art::InputTag> TrackModuleLabel{
        fhicl::Name("TrackModuleLabel"),
          fhicl::Comment("Module that created Tracks"),
          "default"
      };
      fhicl::Atom<art::InputTag> OpHitModuleLabel{
        fhicl::Name("OpHitModuleLabel"),
          fhicl::Comment("Module that created OpHits"),
          "default"
      };
      fhicl::Atom<art::InputTag> FlashModuleLabel{
        fhicl::Name("FlashModuleLabel"),
          fhicl::Comment("Module that created Flashes"),
          "default"
      };
      fhicl::Atom<art::InputTag> OpHitTriggerLabel{
        fhicl::Name("OpHitTriggerLabel"),
          fhicl::Comment("Module that triggered OpHits"),
          "default"
      };
      fhicl::Atom<art::InputTag> FlashTriggerLabel{
        fhicl::Name("FlashTriggerLabel"),
          fhicl::Comment("Module that triggered Flashes"),
          "default"
      };
      fhicl::Atom<art::InputTag> SpacePointLabel{
        fhicl::Name("SpacePointLabel"),
          fhicl::Comment("SpacePoint Solver Module Label"),
          "default"
      };
      fhicl::Atom<art::InputTag> HitModuleLabel{
        fhicl::Name("HitModuleLabel"),
          fhicl::Comment("HitFinder Module Label"),
          "default"
      };
      fhicl::Atom<art::InputTag> ShowerModuleLabel{
        fhicl::Name("ShowerModuleLabel"),
          fhicl::Comment("Shower Finder Module Label"),
          "default"
      };

    };


    StopMuonTracker(fhicl::ParameterSet const& p);
    StopMuonTracker(art::EDAnalyzer::Table<StopMuonTrackerConfig> const& config);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    StopMuonTracker(StopMuonTracker const&) = delete;
    StopMuonTracker(StopMuonTracker&&) = delete;
    StopMuonTracker& operator=(StopMuonTracker const&) = delete;
    StopMuonTracker& operator=(StopMuonTracker&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void beginRun(art::Run const& r) override;
    void endJob() override;
    void endRun(art::Run const& r) override;

  private:

    bool dbg=true;

    art::InputTag fPFPModuleLabel;
    art::InputTag fTrackModuleLabel;
    art::InputTag fOpHitModuleLabel;
    art::InputTag fFlashModuleLabel;
    art::InputTag fOpHitTriggerLabel;
    art::InputTag fFlashTriggerLabel;
    art::InputTag fSpacePointSolverModuleLabel;
    art::InputTag fHitModuleLabel;
    art::InputTag fShowerModuleLabel;

    double trk_len_min_cut = 75.0;
    // Declare member data here.


    ana::smts priv_TrackBranchBuffer;
    Int_t priv_skipped_events_nf_Buffer         =0;
    Int_t priv_skipped_events_oh_Buffer         =0;
    Int_t priv_skipped_events_nt_Buffer         =0;
    Int_t priv_track_has_pfp_Buffer             =0;
    Int_t priv_track_has_t0_Buffer              =0;
    Int_t priv_track_in_bounds_Buffer           =0;
    Int_t priv_track_crosses_cathode_Buffer     =0;
    Int_t priv_track_passes_fiducial_cut_Buffer =0;
    Int_t priv_track_not_broken_Buffer          =0;
    Int_t priv_track_not_short_Buffer           =0;
    Int_t priv_track_has_hits_Buffer            =0;
    Int_t priv_track_has_opflash_Buffer         =0;
    Int_t priv_track_has_no_pfp_Buffer              =0;
    Int_t priv_track_has_no_t0_Buffer               =0;
    Int_t priv_track_not_in_bounds_Buffer           =0;
    Int_t priv_track_no_crosses_cathode_Buffer      =0;
    Int_t priv_track_fails_fiducial_cut_Buffer      =0;
    Int_t priv_track_is_broken_Buffer               =0;
    Int_t priv_track_short_Buffer                   =0;
    Int_t priv_track_has_no_hits_Buffer             =0;
    Int_t priv_track_has_no_ohits_Buffer            =0;
    Int_t priv_track_has_no_oflash_Buffer           =0;
    Int_t priv_track_has_ohits_Buffer               =0;
    TTree* priv_TrackTree=0;
    TTree* priv_StatTree=0;
    TBranch* priv_TrackBranch=0;
    TBranch* priv_StatBranch1=0;
    TBranch* priv_StatBranch2=0;
    TBranch* priv_StatBranch3=0;
    TBranch* priv_StatBranch4=0;
    TBranch* priv_StatBranch5=0;
    TBranch* priv_StatBranch6=0;
    TBranch* priv_StatBranch7=0;
    TBranch* priv_StatBranch8=0;
    TBranch* priv_StatBranch9=0;
    TBranch* priv_StatBranch10=0;
    TBranch* priv_StatBranch11=0;
    TBranch* priv_StatBranch12=0;
    TBranch* priv_StatBranch13=0;
    TBranch* priv_StatBranch14=0;
    TBranch* priv_StatBranch15=0;
    TBranch* priv_StatBranch16=0;
    TBranch* priv_StatBranch17=0;
    TBranch* priv_StatBranch18=0;
    TBranch* priv_StatBranch19=0;
    TBranch* priv_StatBranch20=0;
    TBranch* priv_StatBranch21=0;

};

#endif /* ___StopMuonTracker_h */
