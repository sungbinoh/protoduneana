///////////////////////////////////////////////////////////////////////////////
//
// Author: JStock
// \brief: a short library for different filters and functions related to
// Muon Track Finding
//
//
///////////////////////////////////////////////////////////////////////////////


#include <TVector3.h>

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "dune/DuneObj/SMTracks.h"




///////////////////////////////////////////////////////////////////////////////
//
// Author: JStock
// \brief: Stealing methods from Aleena's Michel Analysis module for use
// in my own work with stop muons.
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
////////////////////Stop Muon Library Headers and Classes//////////////////////
///////////////////////////////////////////////////////////////////////////////


namespace MFL{
  struct MFL_Conf
  {
    public:
      int APAnegbound1() const {return apabound[0];}
      int APAposbound1() const {return apabound[1];}
      int APAnegbound2() const {return apabound[2];}
      int APAposbound2() const {return apabound[3];}
      int apa_bound_sep{10};
      int kMaxCh{90000};
      int kMaxHits{10000};
      double activeBounds_eff[6] = {-359.4155, 359.4155, 0., 607.49875, -0.49375, 695.28625};
      double fThicknessStartVolume = 30;

      double y_end_cut = 80.0;
      double z_end_min_cut = 80.0;
      double z_end_max_cut = 610.0;

      double hitdist=10.0;
      double oltl=10.0;

      int fiducialBounds[6] = {
        int(activeBounds_eff[0])+50,
        int(activeBounds_eff[1])-50,
        int(activeBounds_eff[2])+80,
        int(activeBounds_eff[3])-50,
        int(activeBounds_eff[4])+80,
        int(activeBounds_eff[5])-85
      };
      int fiducialBounds1[6] = {
        int(activeBounds_eff[0])+50,
        int(activeBounds_eff[1])-50,
        int(activeBounds_eff[2])+50,
        int(activeBounds_eff[3])-50,
        int(activeBounds_eff[4])+50,
        int(activeBounds_eff[5])-50
      };
      int apabound[4] =
      {
        int(226-apa_bound_sep),
        int(236+apa_bound_sep),
        int(456-apa_bound_sep),
        int(472+apa_bound_sep)
      };
      double OpHitWindowMin = 50000.0;
      double OpHitWindowMax = 50000.0;
      double OpFlashWindowMin = 50000.0;
      double OpFlashWindowMax = 50000.0;
      double minhitpeakcut = 20.0;
      double maxhitpeakcut = 58000.0;
  };

  const MFL_Conf mfl_conf;

  bool IsPointInBounds(TVector3 const &P);

  bool IsPointInFV(TVector3 const & t);

  void SwitchDirection(TVector3 &ps_s, TVector3 &ps_e, TVector3 &dr_s, TVector3 dr_e);

  void ForceDirection(TVector3 &ps_s, TVector3 &ps_e, TVector3 &dr_s, TVector3 dr_e);

  bool CathCross(TVector3 const& v1, TVector3 const& v2);

  bool FidCut(TVector3 const &p);

  bool BrokenTrack(art::Ptr<recob::Track> const& ptrack, std::vector<art::Ptr<recob::Track>> const& tracklist );

  double MinTime(std::vector<ana::hit> h);

  std::vector<ana::hit> GetHitsForTrack(art::Ptr<recob::Track> ptrack, art::Event const& e, art::InputTag trackModuleLabel); //HitTimesCut is in GetHitsForTrack. Hits sort in GetHitsForTrack.

  //ToDo I htink I should make a vector of ohits instead of a avector of recob::Hits to pass tothis function. The rationale is to make sure that Times are compatible with the TPC.
  std::vector<ana::ohit> GetOpHitsForTrack(art::Ptr<recob::Track> ptrack, double T00, double ophitTrigT, std::vector<art::Ptr<recob::OpHit> > const& ophitList);
  std::vector<ana::oflash> GetOpFlashesForTrack(art::Ptr<recob::Track> ptrack, double T00, double opflashTrigT, std::vector<art::Ptr<recob::OpFlash> > const& opflashList);

  std::vector<ana::hit> DropBadHits(std::vector<ana::hit> & hits, art::FindManyP<recob::Hit, recob::TrackHitMeta> mp_thm, size_t trk_size);
std::vector<std::vector<double>> HasMichel( std::vector<art::Ptr<recob::Track>> const& tracklist, art::FindManyP<recob::Hit, recob::TrackHitMeta> const& fmthm, TVector3 const& pos, TVector3 const& end, size_t const& trknum, double const& T00, art::FindManyP<recob::SpacePoint> const& fmsp, art::Handle<std::vector<recob::Track>> trackListHandle, std::vector<art::Ptr<recob::Hit>> const& hitlist, std::vector<art::Ptr<recob::Shower>> const& showerlist, art::FindManyP<recob::PFParticle> const& pfp_shwr_assn, art::FindManyP<anab::T0> const& shwr_t0_assn_v,  detinfo::DetectorPropertiesData const& detprop );


  //bool HitTimesCut(art::Ptr<recob::Track> const& ptrack, art::Event const& e, std::string label);



}
