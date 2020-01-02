////////////////////////////////////////////////////////////////////////
// Class:       TensionAnalysis
// Plugin Type: analyzer (art v3_01_02)
// File:        TensionAnalysis_module.cc
//
// Generated at Fri Apr 12 14:25:37 2019 by Adam Lister using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "cetlib/search_path.h"

// larsoft
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"

// root
#include "TFile.h"
#include "TTree.h"

// cpp
#include <memory>

// local
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/SelectionCuts.h"
#include "protoduneana/singlephase/tensionanalysis/algorithms/HelperFunctions.h"
#include "protoduneana/singlephase/tensionanalysis/algorithms/Structs.h"

namespace pdsp {
  class TensionAnalysis;
}


class pdsp::TensionAnalysis : public art::EDAnalyzer {
  public:
    explicit TensionAnalysis(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    TensionAnalysis(TensionAnalysis const&) = delete;
    TensionAnalysis(TensionAnalysis&&) = delete;
    TensionAnalysis& operator=(TensionAnalysis const&) = delete;
    TensionAnalysis& operator=(TensionAnalysis&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;

    // used to resize vectors which go into the root tree back to zero.
    // called once per event
    void resizeVectors();

    // given an analysis APA number, return a struct containing the 
    // build APA number along with the side which faces the cathode
    // and therefore should be used for analysis
    pdsp::APAAndSide getAPAInfoFromAnaAPANumber(int anaAPANumber);

    // given a build APA number, return a struct containing the 
    // analysis APA number along with the side which faces the cathode
    // and therefore should be used for analysis
    pdsp::APAAndSide getAPAInfoFromBuildAPANumber(int buildAPANumber);

    // takes a reconstructed spacepoint associated to a hit on a given channel, and 
    // finds the closest wire segment associated with that DAQ channel number
    // then return a the distance to that segment and the tension associated with that 
    // segment
    pdsp::TensionInformation getTensionInformation(art::Ptr< recob::Hit > thisHit, art::Ptr< recob::SpacePoint > thisSpacePoint, TTree* t, pdsp::APASide sideToUse);

    // takes a reconstructes spacepoint associated to a hit on a given channel,
    // and finds the calorimetry information which is located neares to the
    // space point
    pdsp::CalorimetryInformation getCalorimetryInformation(art::Ptr< recob::SpacePoint > thisSpacePoint, std::vector< art::Ptr<anab::Calorimetry> > calVec);

  private:

    // services
    art::ServiceHandle< art::TFileService > tfs;

    // Declare member data here.
    TTree* anaTree;
    ::util::SelectionCuts _selCuts;
    ::util::HelperFunctions _helperFuncs;
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();

    // fhicl parameters
    std::string fPfpLabel;
    std::string fTrackLabel;
    std::string fAllHitLabel;
    std::string fTrackHitLabel;
    std::string fCaloLabel;

    std::string fPfpTrackAssnLabel;
    std::string fTrackT0AssnLabel;
    std::string fTrackHitAssnLabel;
    std::string fHitSpacePointAssnLabel;

    float fMinimumResidualRange;
    float fMaximumResidualRange;
    int   fWireOffsetU;
    int   fWireOffsetV;

    // variables for tree
    int run;
    int subRun;
    int event;
    int isData;

    std::vector<double>* trackStartX                              = nullptr;
    std::vector<double>* trackStartY                              = nullptr;
    std::vector<double>* trackStartZ                              = nullptr;
    std::vector<double>* trackEndX                                = nullptr;
    std::vector<double>* trackEndY                                = nullptr;
    std::vector<double>* trackEndZ                                = nullptr;
    std::vector<double>* trackLength                              = nullptr;
    std::vector<double>* trackTheta                               = nullptr;
    std::vector<double>* trackPhi                                 = nullptr;
    std::vector<double>* trackAzimuthal                           = nullptr;
    std::vector<double>* trackZenith                              = nullptr;
    std::vector<double>* trackThetaXZ                             = nullptr;
    std::vector<double>* trackThetaYZ                             = nullptr;
    std::vector<double>* trackT0                                  = nullptr;
    std::vector< std::vector<float> >* trackHitRMS                = nullptr;
    std::vector< std::vector<float> >* trackHitPeakTime           = nullptr;
    std::vector< std::vector<float> >* trackHitPeakAmplitude      = nullptr;
    std::vector< std::vector<float> >* trackHitIntegral           = nullptr;
    std::vector< std::vector<int> >*   trackHitChannel            = nullptr;
    std::vector< std::vector<int> >*   trackHitView               = nullptr;
    std::vector< std::vector<float> >* trackHitDistanceToTrackEnd = nullptr;
    std::vector< std::vector<int> >*   trackHitAPAAnalysisNumber  = nullptr;
    std::vector< std::vector<int> >*   trackHitAPABuildNumber     = nullptr;
    std::vector< std::vector<int> >*   trackHitWireNo             = nullptr;
    std::vector< std::vector<float> >* trackHitWireLength         = nullptr;
    std::vector< std::vector<float> >* trackHitSegmentDistance    = nullptr;
    std::vector< std::vector<float> >* trackHitWireTension        = nullptr;
    std::vector< std::vector<float> >* trackHitWireSegmentYStart  = nullptr;
    std::vector< std::vector<float> >* trackHitWireSegmentYEnd    = nullptr;
    std::vector< std::vector<float> >* trackHitWireSegmentZStart  = nullptr;
    std::vector< std::vector<float> >* trackHitWireSegmentZEnd    = nullptr;
    std::vector< std::vector<float> >* trackHitWireSegmentLength  = nullptr;
    std::vector< std::vector<float> >* trackHitWireGeomYStart     = nullptr;
    std::vector< std::vector<float> >* trackHitWireGeomYEnd       = nullptr;
    std::vector< std::vector<float> >* trackHitWireGeomZStart     = nullptr;
    std::vector< std::vector<float> >* trackHitWireGeomZEnd       = nullptr;
    std::vector< std::vector<float> >* trackHitWireGeomLength     = nullptr;
    std::vector< std::vector<float> >* trackHitSPCaloDist         = nullptr;
    std::vector< std::vector<float> >* trackHitCaloEnergyDep      = nullptr;
    std::vector< std::vector<float> >* trackHitCaloChargeDep      = nullptr;
    std::vector< std::vector<float> >* trackHitCaloResidualRange  = nullptr;
    std::vector<float>* allHitRMS                                 = nullptr;
    std::vector<float>* allHitPeakTime                            = nullptr;
    std::vector<float>* allHitPeakAmplitude                       = nullptr;
    std::vector<float>* allHitIntegral                            = nullptr;
    std::vector<int>* allHitChannel                               = nullptr;
    std::vector<int>* allHitView                                  = nullptr;

    // histograms to save
    TTree* geometryTree;
    int channelNumber;
    int channelNumberAssociatedWires;
    std::vector<int> channelAssociatedWiresPlane;
    std::vector<int> channelAssociatedAPABuildNumber;
    std::vector<double> channelAssociatedWiresLength;
    std::vector<double> channelAssociatedWiresStartX;
    std::vector<double> channelAssociatedWiresEndX;
    std::vector<double> channelAssociatedWiresStartY;
    std::vector<double> channelAssociatedWiresEndY;
    std::vector<double> channelAssociatedWiresStartZ;
    std::vector<double> channelAssociatedWiresEndZ;

    // root files to read in
    TFile* tensionsFile;

    // TTrees and maps for individual planes
    TTree* treeXLayerUS001;
    TTree* treeULayerUS001;
    TTree* treeVLayerUS001;
    TTree* treeXLayerUS002;
    TTree* treeULayerUS002;
    TTree* treeVLayerUS002;
    TTree* treeXLayerUS003;
    TTree* treeULayerUS003;
    TTree* treeVLayerUS003;
    TTree* treeXLayerUS004;
    TTree* treeULayerUS004;
    TTree* treeVLayerUS004;
    TTree* treeXLayerUK001;
    TTree* treeULayerUK001;
    TTree* treeVLayerUK001;
    TTree* treeXLayerUK002;
    TTree* treeULayerUK002;
    TTree* treeVLayerUK002;

};


pdsp::TensionAnalysis::TensionAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
{

  MF_LOG_VERBATIM("TensionAnalysis")
    << "-- begin TensionAnalysis::TensionAnalysis";

  fhicl::ParameterSet const pLabels = p.get< fhicl::ParameterSet >("ProducerLabels");
  fhicl::ParameterSet const pCuts   = p.get< fhicl::ParameterSet  >("CutValues");

  fPfpLabel               = pLabels.get<std::string>("PFPLabel");
  fTrackLabel             = pLabels.get<std::string>("TrackLabel");
  fAllHitLabel            = pLabels.get<std::string>("AllHitLabel");
  fTrackHitLabel          = pLabels.get<std::string>("TrackHitLabel");
  fCaloLabel              = pLabels.get<std::string>("CalorimetryLabel");
  fPfpTrackAssnLabel      = pLabels.get<std::string>("PFPTrackAssnLabel");
  fTrackT0AssnLabel       = pLabels.get<std::string>("TrackT0AssnLabel");
  fTrackHitAssnLabel      = pLabels.get<std::string>("TrackHitAssnLabel");
  fHitSpacePointAssnLabel = pLabels.get<std::string>("HitSpacePointAssnLabel");

  fMinimumResidualRange   = pCuts.get<float>("MinimumResidualRange");
  fMaximumResidualRange   = pCuts.get<float>("MaximumResidualRange");
  fWireOffsetU            = pCuts.get<int>  ("WireOffsetU");
  fWireOffsetV            = pCuts.get<int>  ("WireOffsetV");

  _selCuts.Configure(pCuts);

  MF_LOG_VERBATIM("TensionAnalysis")
    << "---------- Printing Configuration ------------"
    << "\n -- fPfpLabel               : " << fPfpLabel
    << "\n -- fTrackLabel             : " << fTrackLabel
    << "\n -- fAllHitLabel            : " << fAllHitLabel
    << "\n -- fTrackHitLabel          : " << fTrackHitLabel
    << "\n -- fCaloLabel              : " << fCaloLabel
    << "\n -- fPfpTrackAssnLabel      : " << fPfpTrackAssnLabel
    << "\n -- fTrackT0AssnLabel       : " << fTrackT0AssnLabel
    << "\n -- fTrackHitLabel          : " << fTrackHitLabel
    << "\n -- fHitSpacePointAssnLabel : " << fHitSpacePointAssnLabel
    << "\n -- fMinimumResidualRange   : " << fMinimumResidualRange
    << "\n -- fMaximumResidualRange   : " << fMaximumResidualRange
    << "\n -- fWireOffsetU            : " << fWireOffsetU
    << "\n -- fWireOffsetV            : " << fWireOffsetV
    << "\n----------------------------------------------";

}

void pdsp::TensionAnalysis::analyze(art::Event const& e)
{

  this->resizeVectors();

  // get auxiliary information
  run    = e.run();
  subRun = e.subRun();
  event  = e.event();
  isData = e.isRealData();

  MF_LOG_VERBATIM("TensionAnalysis")
    << "Processing event " << run << "." << subRun << "." << event; 

  // get handles to information we're interested in
  art::Handle< std::vector< recob::PFParticle > > pfpHandle;
  e.getByLabel(fPfpLabel, pfpHandle);
  std::vector< art::Ptr< recob::PFParticle > > pfpPtrVector;
  art::fill_ptr_vector(pfpPtrVector, pfpHandle);

  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector< art::Ptr< recob::Track > > trackPtrVector;
  art::fill_ptr_vector(trackPtrVector, trackHandle);

  // need to have two hit collections
  // -- allHits is for all hits in the event
  // -- trackHits is just those associated with selected tracks
  //    (allows access to association)

  art::Handle< std::vector< recob::Hit > > allHitHandle;
  e.getByLabel(fAllHitLabel, allHitHandle);
  std::vector< art::Ptr< recob::Hit > > allHitPtrVector;
  art::fill_ptr_vector(allHitPtrVector, allHitHandle);

  art::Handle< std::vector< recob::Hit > > trackHitHandle;
  e.getByLabel(fTrackHitLabel, trackHitHandle);
  std::vector< art::Ptr< recob::Hit > > trackHitPtrVector;
  art::fill_ptr_vector(trackHitPtrVector, trackHitHandle);

  MF_LOG_DEBUG("pdsp::TensionAnalysis")
  << "getting calorimetry information...";
  art::Handle< std::vector< anab::Calorimetry > > caloHandle;
  e.getByLabel(fCaloLabel, caloHandle);
  std::vector< art::Ptr< anab::Calorimetry > > caloPtrVector;
  art::fill_ptr_vector(caloPtrVector, caloHandle);
  MF_LOG_DEBUG("pdsp::TensionAnalysis")
  << "got calorimetry information...";

  // associations
  art::FindManyP< recob::Track >      tracksFromPfps(pfpHandle, e, fPfpTrackAssnLabel);
  art::FindManyP< anab::T0 >          t0sFromTracks(trackHandle, e, fTrackT0AssnLabel);
  art::FindManyP< recob::Hit >        hitsFromTracks(trackHandle, e, fTrackHitAssnLabel);
  art::FindManyP< recob::SpacePoint > spacePointFromHits(trackHitHandle, e, fHitSpacePointAssnLabel);

  // loop hits for dead channel analysis
  for (size_t i = 0; i < allHitPtrVector.size(); i++){

    art::Ptr< recob::Hit > thisHit = allHitPtrVector.at(i);
    allHitRMS          ->push_back(thisHit->RMS());
    allHitPeakTime     ->push_back(thisHit->PeakTime());
    allHitPeakAmplitude->push_back(thisHit->PeakAmplitude());
    allHitIntegral     ->push_back(thisHit->Integral());
    allHitChannel      ->push_back(thisHit->Channel());
    allHitView         ->push_back(thisHit->View());

  }

  // loop pfp for tension analysis
  for ( size_t i = 0; i < pfpPtrVector.size(); i++){

    art::Ptr< recob::PFParticle > thisPfp = pfpPtrVector.at(i);

    if ((tracksFromPfps.at(thisPfp.key())).size() == 0) continue;

    // find particles which are identified as track-like
    if (thisPfp->PdgCode() == 13){

      // get associated tracks
      art::Ptr< recob::Track > thisTrack = (tracksFromPfps.at(thisPfp.key())).at(0);

      trackStartX   ->push_back(thisTrack->Start().X());
      trackStartY   ->push_back(thisTrack->Start().Y());
      trackStartZ   ->push_back(thisTrack->Start().Z());
      trackEndX     ->push_back(thisTrack->End().X());
      trackEndY     ->push_back(thisTrack->End().Y());
      trackEndZ     ->push_back(thisTrack->End().Z());
      trackLength   ->push_back(thisTrack->Length());
      trackTheta    ->push_back(thisTrack->Theta());
      trackPhi      ->push_back(thisTrack->Phi());
      trackAzimuthal->push_back(thisTrack->AzimuthAngle());
      trackZenith   ->push_back(thisTrack->ZenithAngle());

      // this is just reusing code from the stopping muon selection,
      // doesn't actually perform cuts
      std::vector< pdsp::CutResult > isPassThetaXZ = _selCuts.IsPassesThetaXZSelection(thisTrack);
      std::vector< pdsp::CutResult > isPassThetaYZ = _selCuts.IsPassesThetaYZSelection(thisTrack);
      trackThetaXZ->push_back((isPassThetaXZ.at(2)).value);
      trackThetaYZ->push_back((isPassThetaYZ.at(2)).value);

      // get hits associated to this track
      std::vector<art::Ptr< recob::Hit > > theseHits = hitsFromTracks.at(thisTrack.key());

      std::vector<float> tHitRMS;
      std::vector<float> tHitPeakTime;
      std::vector<float> tHitPeakAmplitude;
      std::vector<float> tHitIntegral;
      std::vector<int>   tHitChannel;
      std::vector<int>   tHitView;
      std::vector<float> tHitDistanceToTrackEnd;
      std::vector<int>   tHitAPAAnalysisNumber;
      std::vector<int>   tHitAPABuildNumber;
      std::vector<int>   tHitWireNo;
      std::vector<float> tHitWireLength;
      std::vector<float> tHitSegmentDistance;
      std::vector<float> tHitWireTension;
      std::vector<float> tHitWireSegmentYStart;
      std::vector<float> tHitWireSegmentYEnd;
      std::vector<float> tHitWireSegmentZStart;
      std::vector<float> tHitWireSegmentZEnd;
      std::vector<float> tHitWireSegmentLength;
      std::vector<float> tHitWireGeomYStart;
      std::vector<float> tHitWireGeomYEnd;
      std::vector<float> tHitWireGeomZStart;
      std::vector<float> tHitWireGeomZEnd;
      std::vector<float> tHitWireGeomLength;
      std::vector<float> tHitSPCaloDist;
      std::vector<float> tHitCaloEnergyDep;
      std::vector<float> tHitCaloChargeDep;
      std::vector<float> tHitCaloResidualRange;

      MF_LOG_DEBUG("pdsp::TensionAnalysis::analyze")
        << "-- found " << theseHits.size() << " associated hits. looping.";

      for (size_t iHit = 0; iHit < theseHits.size(); iHit++){
        art::Ptr< recob::Hit > thisHit = theseHits.at(iHit);

        if (spacePointFromHits.at(thisHit.key()).size() == 0){
          MF_LOG_DEBUG("pdsp::TensionAnalysis::analyze") 
            << "-- -- no associated spacepoints, skipping hit.";
          continue;
        }
        art::Ptr< recob::SpacePoint > thisSpacePoint = spacePointFromHits.at(thisHit.key()).at(0);

        double distanceToTrackEnd = 
          _helperFuncs.CalculateDistance({thisSpacePoint->XYZ()[0], thisSpacePoint->XYZ()[1], thisSpacePoint->XYZ()[2]},
              {thisTrack->End().X(), thisTrack->End().Y(), thisTrack->End().Z()});

        // locate hits inside the MIP-region
        if (distanceToTrackEnd > fMinimumResidualRange && distanceToTrackEnd < fMaximumResidualRange){
          tHitRMS          .push_back(thisHit->RMS());
          tHitPeakTime     .push_back(thisHit->PeakTime());
          tHitPeakAmplitude.push_back(thisHit->PeakAmplitude());
          tHitIntegral     .push_back(thisHit->Integral());
          tHitChannel      .push_back(thisHit->Channel());
          tHitView         .push_back(thisHit->View());
          tHitDistanceToTrackEnd.push_back(distanceToTrackEnd);

          // for this hit, fill an pdsp::APAAndSide object, which stores information
          // about the APA, and which side of the APA is used for analysis
          int hitAPAAnalysisNumber = thisHit->WireID().asTPCID().TPC;
          pdsp::APAAndSide thisAPAAndSide = this->getAPAInfoFromAnaAPANumber(hitAPAAnalysisNumber);

          tHitAPAAnalysisNumber.push_back(thisAPAAndSide.APAAnalysisNumber);
          tHitAPABuildNumber   .push_back(thisAPAAndSide.APABuildNumber);

          // go through the geometry in order to find the the WireID, along with
          // its length
          geo::WireID const& thisWire = thisHit->WireID();
          int hitWireNo = thisWire.Wire;
          tHitWireNo.push_back(hitWireNo);

          geo::WireGeo const& thisWireGeo = geom->Wire(thisWire);
          tHitWireLength.push_back(thisWireGeo.Length());

          MF_LOG_DEBUG("TensionAnalysis::analyze") 
            << "hit has TPC number: " << hitAPAAnalysisNumber 
            << "  which corresponds to build number : " << thisAPAAndSide.APABuildNumber;


          // now, dependent on what the APA and side is, 
          // then we need to access a different tree to pull out the tension information
          // there are definitely smarter ways we can deal with this but for now, let's be explicit

          pdsp::TensionInformation thisTensionInformation;
          
          if (thisAPAAndSide.APABuildNumber == 1){
            MF_LOG_DEBUG("TensionAnalysis::analyze") << "using Build APA 1";
            if (thisHit->View() == 2)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeXLayerUS001, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 1)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeVLayerUS001, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 0)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeULayerUS001, thisAPAAndSide.sideToUse);
          }
          if (thisAPAAndSide.APABuildNumber == 2){
            MF_LOG_DEBUG("TensionAnalysis::analyze") << "using Build APA 2";
            if (thisHit->View() == 2)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeXLayerUS002, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 1)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeVLayerUS002, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 0)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeULayerUS002, thisAPAAndSide.sideToUse);
          }
          if (thisAPAAndSide.APABuildNumber == 3){
            MF_LOG_DEBUG("TensionAnalysis::analyze") << "using Build APA 3";
            if (thisHit->View() == 2)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeXLayerUS003, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 1)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeVLayerUS003, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 0)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeULayerUS003, thisAPAAndSide.sideToUse);
          }
          if (thisAPAAndSide.APABuildNumber == 4){
            MF_LOG_DEBUG("TensionAnalysis::analyze") << "using Build APA 4";
            if (thisHit->View() == 2)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeXLayerUS004, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 1)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeVLayerUS004, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 0)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeULayerUS004, thisAPAAndSide.sideToUse);
          }
          if (thisAPAAndSide.APABuildNumber == 5){
            MF_LOG_DEBUG("TensionAnalysis::analyze") << "using Build APA 5";
            if (thisHit->View() == 2)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeXLayerUK001, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 1)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeVLayerUK001, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 0)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeULayerUK001, thisAPAAndSide.sideToUse);
          }
          if (thisAPAAndSide.APABuildNumber == 6){
            MF_LOG_DEBUG("TensionAnalysis::analyze") << "using Build APA 6";
            if (thisHit->View() == 2)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeXLayerUK002, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 1)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeVLayerUK002, thisAPAAndSide.sideToUse);
            if (thisHit->View() == 0)
              thisTensionInformation = this->getTensionInformation(thisHit, thisSpacePoint, treeULayerUK002, thisAPAAndSide.sideToUse);
          }

          // now pull out information fro the tension information object to fill the tree
          tHitSegmentDistance   .push_back(thisTensionInformation.distSPGeometry     );
          tHitWireTension       .push_back(thisTensionInformation.tension            );
          tHitWireGeomYStart    .push_back(thisTensionInformation.wireGeometryYStart );
          tHitWireGeomYEnd      .push_back(thisTensionInformation.wireGeometryYEnd   );
          tHitWireGeomZStart    .push_back(thisTensionInformation.wireGeometryZStart );
          tHitWireGeomZEnd      .push_back(thisTensionInformation.wireGeometryZEnd   );
          tHitWireGeomLength    .push_back(thisTensionInformation.wireGeometryLength );
          tHitWireSegmentYStart .push_back(thisTensionInformation.wireSegmentYStart  );
          tHitWireSegmentYEnd   .push_back(thisTensionInformation.wireSegmentYEnd    );
          tHitWireSegmentZStart .push_back(thisTensionInformation.wireSegmentZStart  );
          tHitWireSegmentZEnd   .push_back(thisTensionInformation.wireSegmentZEnd    );
          tHitWireSegmentLength .push_back(thisTensionInformation.wireSegmentLength  );

          // get calorimetry information for this spacepoint
          pdsp::CalorimetryInformation thisCalorimetryInfo = this->getCalorimetryInformation(thisSpacePoint, caloPtrVector);
          tHitSPCaloDist       .push_back(thisCalorimetryInfo.distanceFromSPToCalo);
          tHitCaloEnergyDep    .push_back(thisCalorimetryInfo.energyDeposition);
          tHitCaloChargeDep    .push_back(thisCalorimetryInfo.chargeDeposition);
          tHitCaloResidualRange.push_back(thisCalorimetryInfo.residualRange);

        }
        else {
          MF_LOG_DEBUG("pdsp::TensionAnalysis::analyze")
            << "-- -- hit outside of MIP region, skipping";
          tHitSegmentDistance   .push_back(-1);
          tHitWireTension       .push_back(-1);
          tHitWireGeomYStart    .push_back(-1);
          tHitWireGeomYEnd      .push_back(-1);
          tHitWireGeomZStart    .push_back(-1);
          tHitWireGeomZEnd      .push_back(-1);
          tHitWireGeomLength    .push_back(-1);
          tHitWireSegmentYStart .push_back(-1);
          tHitWireSegmentYEnd   .push_back(-1);
          tHitWireSegmentZStart .push_back(-1);
          tHitWireSegmentZEnd   .push_back(-1);
          tHitWireSegmentLength .push_back(-1);
          tHitSPCaloDist        .push_back(-1);
          tHitCaloEnergyDep     .push_back(-1);
          tHitCaloChargeDep     .push_back(-1);
          tHitCaloResidualRange .push_back(-1);

        }

      }

      trackHitRMS               -> push_back(tHitRMS);
      trackHitPeakTime          -> push_back(tHitPeakTime);
      trackHitPeakAmplitude     -> push_back(tHitPeakAmplitude);
      trackHitIntegral          -> push_back(tHitIntegral);
      trackHitChannel           -> push_back(tHitChannel);
      trackHitView              -> push_back(tHitView);
      trackHitDistanceToTrackEnd-> push_back(tHitDistanceToTrackEnd);
      trackHitAPAAnalysisNumber -> push_back(tHitAPAAnalysisNumber);
      trackHitAPABuildNumber    -> push_back(tHitAPABuildNumber);
      trackHitWireNo            -> push_back(tHitWireNo);
      trackHitWireLength        -> push_back(tHitWireLength);
      trackHitSegmentDistance   -> push_back(tHitSegmentDistance);
      trackHitWireTension       -> push_back(tHitWireTension);
      trackHitWireSegmentYStart -> push_back(tHitWireSegmentYStart);
      trackHitWireSegmentYEnd   -> push_back(tHitWireSegmentYEnd);
      trackHitWireSegmentZStart -> push_back(tHitWireSegmentZStart);
      trackHitWireSegmentZEnd   -> push_back(tHitWireSegmentZEnd);
      trackHitWireSegmentLength -> push_back(tHitWireSegmentLength);
      trackHitWireGeomYStart    -> push_back(tHitWireGeomYStart);
      trackHitWireGeomYEnd      -> push_back(tHitWireGeomYEnd);
      trackHitWireGeomZStart    -> push_back(tHitWireGeomZStart);
      trackHitWireGeomZEnd      -> push_back(tHitWireGeomZEnd);
      trackHitWireGeomLength    -> push_back(tHitWireGeomLength);
      trackHitSPCaloDist        -> push_back(tHitSPCaloDist);
      trackHitCaloEnergyDep     -> push_back(tHitCaloEnergyDep);
      trackHitCaloChargeDep     -> push_back(tHitCaloChargeDep);
      trackHitCaloResidualRange -> push_back(tHitCaloResidualRange);
      

      art::Ptr< anab::T0 > thisT0 = t0sFromTracks.at(thisTrack.key()).at(0);

      trackT0->push_back(thisT0->Time());

    }

  }

  anaTree->Fill();
}

void pdsp::TensionAnalysis::beginJob()
{
  MF_LOG_DEBUG("TensionAnalysis")
    << "-- begin TensionAnalysis::beginJob";

  // Implementation of optional member functon here.
  anaTree = tfs->make<TTree>("analysis_tree"  , "analysis tree");
  anaTree->Branch("run"                       , &run);
  anaTree->Branch("subRun"                    , &subRun);
  anaTree->Branch("event"                     , &event);
  anaTree->Branch("trackStartX"               , "std::vector< double >"             , &trackStartX);
  anaTree->Branch("trackStartY"               , "std::vector< double >"             , &trackStartY);
  anaTree->Branch("trackStartZ"               , "std::vector< double >"             , &trackStartZ);
  anaTree->Branch("trackEndX"                 , "std::vector< double >"             , &trackEndX);
  anaTree->Branch("trackEndY"                 , "std::vector< double >"             , &trackEndY);
  anaTree->Branch("trackEndZ"                 , "std::vector< double >"             , &trackEndZ);
  anaTree->Branch("trackLength"               , "std::vector< double >"             , &trackLength);
  anaTree->Branch("trackTheta"                , "std::vector< double >"             , &trackTheta);
  anaTree->Branch("trackPhi"                  , "std::vector< double >"             , &trackPhi);
  anaTree->Branch("trackAzimuthal"            , "std::vector< double >"             , &trackAzimuthal);
  anaTree->Branch("trackZenith"               , "std::vector< double >"             , &trackZenith);
  anaTree->Branch("trackThetaXZ"              , "std::vector< double >"             , &trackThetaXZ);
  anaTree->Branch("trackThetaYZ"              , "std::vector< double >"             , &trackThetaYZ);
  anaTree->Branch("trackT0"                   , "std::vector< double >"             , &trackT0);
  anaTree->Branch("trackHitRMS"               , "std::vector< std::vector<float> >" , &trackHitRMS);
  anaTree->Branch("trackHitPeakTime"          , "std::vector< std::vector<float> >" , &trackHitPeakTime);
  anaTree->Branch("trackHitPeakAmplitude"     , "std::vector< std::vector<float> >" , &trackHitPeakAmplitude);
  anaTree->Branch("trackHitIntegral"          , "std::vector< std::vector<float> >" , &trackHitIntegral);
  anaTree->Branch("trackHitChannel"           , "std::vector< std::vector<int> >"   , &trackHitChannel);
  anaTree->Branch("trackHitView"              , "std::vector< std::vector<int> >"   , &trackHitView);
  anaTree->Branch("trackHitDistanceToTrackEnd", "std::vector< std::vector<float> >" , &trackHitDistanceToTrackEnd);
  anaTree->Branch("trackHitAPAAnalysisNumber" , "std::vector< std::vector<int> >"   , &trackHitAPAAnalysisNumber);
  anaTree->Branch("trackHitAPABuildNumber"    , "std::vector< std::vector<int> >"   , &trackHitAPABuildNumber);
  anaTree->Branch("trackHitWireNo"            , "std::vector< std::vector<int> >"   , &trackHitWireNo);
  anaTree->Branch("trackHitWireLength"        , "std::vector< std::vector<float> >" , &trackHitWireLength);
  anaTree->Branch("trackHitSegmentDistance"   , "std::vector< std::vector<float> >" , &trackHitSegmentDistance);
  anaTree->Branch("trackHitWireTension"       , "std::vector< std::vector<float> >" , &trackHitWireTension);
  anaTree->Branch("trackHitWireSegmentYStart" , "std::vector< std::vector<float> >" , &trackHitWireSegmentYStart);
  anaTree->Branch("trackHitWireSegmentYEnd"   , "std::vector< std::vector<float> >" , &trackHitWireSegmentYEnd);
  anaTree->Branch("trackHitWireSegmentZStart" , "std::vector< std::vector<float> >" , &trackHitWireSegmentZStart);
  anaTree->Branch("trackHitWireSegmentZEnd"   , "std::vector< std::vector<float> >" , &trackHitWireSegmentZEnd);
  anaTree->Branch("trackHitWireSegmentLength" , "std::vector< std::vector<float> >" , &trackHitWireSegmentLength);
  anaTree->Branch("trackHitWireGeomYStart"    , "std::vector< std::vector<float> >" , &trackHitWireGeomYStart);
  anaTree->Branch("trackHitWireGeomYEnd"      , "std::vector< std::vector<float> >" , &trackHitWireGeomYEnd);
  anaTree->Branch("trackHitWireGeomZStart"    , "std::vector< std::vector<float> >" , &trackHitWireGeomZStart);
  anaTree->Branch("trackHitWireGeomZEnd"      , "std::vector< std::vector<float> >" , &trackHitWireGeomZEnd);
  anaTree->Branch("trackHitWireGeomLength"    , "std::vector< std::vector<float> >" , &trackHitWireGeomLength);
  anaTree->Branch("trackHitSPCaloDist"        , "std::vector< std::vector<float> >" , &trackHitSPCaloDist);
  anaTree->Branch("trackHitCaloEnergyDep"     , "std::vector< std::vector<float> >" , &trackHitCaloEnergyDep);
  anaTree->Branch("trackHitCaloChargeDep"     , "std::vector< std::vector<float> >" , &trackHitCaloChargeDep);
  anaTree->Branch("trackHitCaloResidualRange" , "std::vector< std::vector<float> >" , &trackHitCaloResidualRange);
  anaTree->Branch("allHitRMS"                 , "std::vector<float>"                , &allHitRMS);
  anaTree->Branch("allHitPeakTime"            , "std::vector<float>"                , &allHitPeakTime);
  anaTree->Branch("allHitPeakAmplitude"       , "std::vector<float>"                , &allHitPeakAmplitude);
  anaTree->Branch("allHitIntegral"            , "std::vector<float>"                , &allHitIntegral);
  anaTree->Branch("allHitChannel"             , "std::vector<int>"                  , &allHitChannel);
  anaTree->Branch("allHitView"                , "std::vector<int>"                  , &allHitView);

  geometryTree = tfs->make<TTree>("geometry_tree"        , "geometry tree");
  geometryTree->Branch("channelNumber"                   , &channelNumber);
  geometryTree->Branch("channelNumberAssociatedWires"    , &channelNumberAssociatedWires);
  geometryTree->Branch("channelAssociatedWiresPlane"     , &channelAssociatedWiresPlane);
  geometryTree->Branch("channelAssociatedAPABuildNumber" , &channelAssociatedAPABuildNumber);
  geometryTree->Branch("channelAssociatedWiresLength"    , "std::vector<double>"              , &channelAssociatedWiresLength);
  geometryTree->Branch("channelAssociatedWiresStartX"    , "std::vector<double>"              , &channelAssociatedWiresStartX);
  geometryTree->Branch("channelAssociatedWiresEndX"      , "std::vector<double>"              , &channelAssociatedWiresEndX);
  geometryTree->Branch("channelAssociatedWiresStartY"    , "std::vector<double>"              , &channelAssociatedWiresStartY);
  geometryTree->Branch("channelAssociatedWiresEndY"      , "std::vector<double>"              , &channelAssociatedWiresEndY);
  geometryTree->Branch("channelAssociatedWiresStartZ"    , "std::vector<double>"              , &channelAssociatedWiresStartZ);
  geometryTree->Branch("channelAssociatedWiresEndZ"      , "std::vector<double>"              , &channelAssociatedWiresEndZ);

  // read in ROOT trees containing wire information
  std::string tensionPath;
  cet::search_path sp("FW_SEARCH_PATH");
  if(!sp.find_file("tensionanalysis/data/tension_measurements.root", tensionPath)){
    throw cet::exception("FileError")
      << "Cannot find tension_measurements.root file "
      << " bail ungracefully\n\n"
      << __FILE__ << ":" << __LINE__;
  }
  tensionsFile = new TFile(tensionPath.c_str(), "read");

  treeXLayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_XLAYER");
  treeULayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_ULAYER");
  treeVLayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_VLAYER");
  treeXLayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_XLAYER");
  treeULayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_ULAYER");
  treeVLayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_VLAYER");
  treeXLayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_XLAYER");
  treeULayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_ULAYER");
  treeVLayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_VLAYER");
  treeXLayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_XLAYER");
  treeULayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_ULAYER");
  treeVLayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_VLAYER");
  treeXLayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_XLAYER");
  treeULayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_ULAYER");
  treeVLayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_VLAYER");
  treeXLayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_XLAYER");
  treeULayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_ULAYER");
  treeVLayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_VLAYER");

  // loop all channels and get the number of assoicated wires, along with their lengths
  for (int i_chan = 0; i_chan < 15360; i_chan++){
    channelAssociatedWiresLength.resize(0);
    channelAssociatedAPABuildNumber.resize(0);
    channelAssociatedWiresPlane.resize(0);
    channelAssociatedWiresStartX.resize(0);
    channelAssociatedWiresStartY.resize(0);
    channelAssociatedWiresStartZ.resize(0);
    channelAssociatedWiresEndX.resize(0);
    channelAssociatedWiresEndY.resize(0);
    channelAssociatedWiresEndZ.resize(0);
    std::vector<geo::WireID> wireIDs = geom->ChannelToWire(i_chan);

    channelNumber = i_chan;
    channelNumberAssociatedWires = (int)wireIDs.size();

    for (size_t i_wire = 0; i_wire < wireIDs.size(); i_wire++){
      geo::WireID thisWireID = wireIDs.at(i_wire);
      geo::WireGeo const& thisWireGeo = geom->Wire(thisWireID);
      channelAssociatedWiresLength.push_back(thisWireGeo.Length());
      channelAssociatedWiresPlane.push_back(thisWireID.Plane);
      channelAssociatedWiresStartX.push_back(thisWireGeo.GetStart()[0]);
      channelAssociatedWiresEndX.push_back(thisWireGeo.GetEnd()[0]);
      channelAssociatedWiresStartY.push_back(thisWireGeo.GetStart()[1]);
      channelAssociatedWiresEndY.push_back(thisWireGeo.GetEnd()[1]);
      channelAssociatedWiresStartZ.push_back(thisWireGeo.GetStart()[2]);
      channelAssociatedWiresEndZ.push_back(thisWireGeo.GetEnd()[2]);


      // get TPC ID information
      pdsp::APAAndSide thisAPAAndSide = this->getAPAInfoFromAnaAPANumber(thisWireID.TPC);
      channelAssociatedAPABuildNumber.push_back(thisAPAAndSide.APABuildNumber);

    }

    geometryTree->Fill();
  }

  MF_LOG_DEBUG("TensionAnalysis")
    << "-- end TensionAnalysis::beginJob";

}

void pdsp::TensionAnalysis::resizeVectors()
{

  MF_LOG_DEBUG("pdsp::TensionAnalysis::resizeVectors")
    << "begin resizeVectors";

  // resize variables
  trackStartX               -> resize(0);
  trackStartY               -> resize(0);
  trackStartZ               -> resize(0);
  trackEndX                 -> resize(0);
  trackEndY                 -> resize(0);
  trackEndZ                 -> resize(0);
  trackLength               -> resize(0);
  trackTheta                -> resize(0);
  trackPhi                  -> resize(0);
  trackAzimuthal            -> resize(0);
  trackZenith               -> resize(0);
  trackThetaXZ              -> resize(0);
  trackThetaYZ              -> resize(0);
  trackT0                   -> resize(0);
  trackHitRMS               -> resize(0);
  trackHitPeakTime          -> resize(0);
  trackHitPeakAmplitude     -> resize(0);
  trackHitIntegral          -> resize(0);
  trackHitChannel           -> resize(0);
  trackHitView              -> resize(0);
  trackHitDistanceToTrackEnd-> resize(0);
  trackHitAPAAnalysisNumber -> resize(0);
  trackHitAPABuildNumber    -> resize(0);
  trackHitWireNo            -> resize(0);
  trackHitWireLength        -> resize(0);
  trackHitSegmentDistance   -> resize(0);
  trackHitWireTension       -> resize(0);
  trackHitWireSegmentYStart -> resize(0);
  trackHitWireSegmentYEnd   -> resize(0);
  trackHitWireSegmentZStart -> resize(0);
  trackHitWireSegmentZEnd   -> resize(0);
  trackHitWireSegmentLength -> resize(0);
  trackHitWireGeomYStart    -> resize(0);
  trackHitWireGeomYEnd      -> resize(0);
  trackHitWireGeomZStart    -> resize(0);
  trackHitWireGeomZEnd      -> resize(0);
  trackHitWireGeomLength    -> resize(0);
  trackHitSPCaloDist        -> resize(0);
  trackHitCaloEnergyDep     -> resize(0);
  trackHitCaloChargeDep     -> resize(0);
  trackHitCaloResidualRange -> resize(0);
  allHitRMS                 -> resize(0);
  allHitPeakTime            -> resize(0);
  allHitPeakAmplitude       -> resize(0);
  allHitIntegral            -> resize(0);
  allHitChannel             -> resize(0);
  allHitView                -> resize(0);

  MF_LOG_DEBUG("pdsp::TensionAnalysis::resizeVectors")
    << "end resizeVectors";

}

pdsp::APAAndSide pdsp::TensionAnalysis::getAPAInfoFromBuildAPANumber(int buildAPANumber)
{

  pdsp::APAAndSide thisAPAAndSide;

  if (buildAPANumber == 1){
    thisAPAAndSide.APABuildNumber         =  1;
    thisAPAAndSide.APAInstallationNumber  =  1;
    thisAPAAndSide.APAAnalysisNumber      =  9;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (buildAPANumber == 2){
    thisAPAAndSide.APABuildNumber         =  2;
    thisAPAAndSide.APAInstallationNumber  =  2;
    thisAPAAndSide.APAAnalysisNumber      =  5;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (buildAPANumber == 3){
    thisAPAAndSide.APABuildNumber         =  3;
    thisAPAAndSide.APAInstallationNumber  =  4;
    thisAPAAndSide.APAAnalysisNumber      =  10;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else if (buildAPANumber == 4){
    thisAPAAndSide.APABuildNumber         =  4;
    thisAPAAndSide.APAInstallationNumber  =  6;
    thisAPAAndSide.APAAnalysisNumber      =  6;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else if (buildAPANumber == 5){
    thisAPAAndSide.APABuildNumber         =  5;
    thisAPAAndSide.APAInstallationNumber  =  3;
    thisAPAAndSide.APAAnalysisNumber      =  1;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (buildAPANumber == 6){
    thisAPAAndSide.APABuildNumber         =  6;
    thisAPAAndSide.APAInstallationNumber  =  5;
    thisAPAAndSide.APAAnalysisNumber      =  2;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else {
    thisAPAAndSide.APABuildNumber         =  -1;
    thisAPAAndSide.APAInstallationNumber  =  -1;
    thisAPAAndSide.APAAnalysisNumber      =  -1;
    thisAPAAndSide.sideToUse              =  pdsp::kUnknown;
  }

  return thisAPAAndSide;
}

pdsp::APAAndSide pdsp::TensionAnalysis::getAPAInfoFromAnaAPANumber(int anaAPANumber)
{

  pdsp::APAAndSide thisAPAAndSide;

  if (anaAPANumber == 1){
    thisAPAAndSide.APABuildNumber         =  5;
    thisAPAAndSide.APAInstallationNumber  =  3;
    thisAPAAndSide.APAAnalysisNumber      =  1;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (anaAPANumber == 2){
    thisAPAAndSide.APABuildNumber         =  6;
    thisAPAAndSide.APAInstallationNumber  =  5;
    thisAPAAndSide.APAAnalysisNumber      =  2;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else if (anaAPANumber == 5){
    thisAPAAndSide.APABuildNumber         =  2;
    thisAPAAndSide.APAInstallationNumber  =  2;
    thisAPAAndSide.APAAnalysisNumber      =  5;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (anaAPANumber == 6){
    thisAPAAndSide.APABuildNumber         =  4;
    thisAPAAndSide.APAInstallationNumber  =  6;
    thisAPAAndSide.APAAnalysisNumber      =  6;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else if (anaAPANumber == 9){
    thisAPAAndSide.APABuildNumber         =  1;
    thisAPAAndSide.APAInstallationNumber  =  1;
    thisAPAAndSide.APAAnalysisNumber      =  9;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (anaAPANumber == 10){
    thisAPAAndSide.APABuildNumber         =  3;
    thisAPAAndSide.APAInstallationNumber  =  4;
    thisAPAAndSide.APAInstallationNumber  =  10;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else {
    thisAPAAndSide.APABuildNumber         =  -1;
    thisAPAAndSide.APAInstallationNumber  =  -1;
    thisAPAAndSide.APAAnalysisNumber      =  -1;
    thisAPAAndSide.sideToUse              =  pdsp::kUnknown;
  }

  return thisAPAAndSide;

}

pdsp::TensionInformation pdsp::TensionAnalysis::getTensionInformation(art::Ptr< recob::Hit > thisHit, art::Ptr< recob::SpacePoint > thisSpacePoint, TTree* t, pdsp::APASide sideToUse){

  MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
    << "beginning search for closest segment";

  pdsp::TensionInformation thisTensionInformation; 

  int   segmentChannel;
  int   segmentNumber;
  float segmentTension;
  float segmentStartY;
  float segmentEndY;
  float segmentStartZ;
  float segmentEndZ;

  int wireOffset = 0;
  if (thisHit->View() == 0)
    wireOffset = fWireOffsetU;
  if (thisHit->View() == 1)
    wireOffset = fWireOffsetV;

  std::string thisSide;
  if (sideToUse == pdsp::kA)
    thisSide = "side_a_";
  else if (sideToUse == pdsp::kB)
    thisSide = "side_b_";
  else throw std::logic_error("pdsp::TensionAnalysis::getTensionInformation --- pdsp::APASide is invalid");

  // numbers from the excel sheets use a different co-ordinate system
  // the APA measurements were taken with the APA in the horizontal orientation
  // meaning that x measures the the longest edge, and y measures the shortest edge
  // 
  // when querying the geometry, z goes along the direction of the shortest edge and
  // y measures the longest edge, so rename variables here.

  t->SetBranchAddress((thisSide+std::string("channel_number")).c_str() , &segmentChannel);
  t->SetBranchAddress((thisSide+std::string("final_tension")).c_str()  , &segmentTension);
  t->SetBranchAddress("segment_number"                                 , &segmentNumber);
  t->SetBranchAddress("x_start"                                        , &segmentStartY);
  t->SetBranchAddress("x_end"                                          , &segmentEndY);
  t->SetBranchAddress("y_start"                                        , &segmentStartZ);
  t->SetBranchAddress("y_end"                                          , &segmentEndZ);

  MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
    << "-- set branch addresses successfully";

  // loop over the tree, find segments corresponding to the channel numbers of the hit
  // and then save the start, end positions of the segment, as well as the tension 
  // of the segment

  std::vector<float> daqChannelWireGeomStartY;
  std::vector<float> daqChannelWireGeomEndY;
  std::vector<float> daqChannelWireGeomStartZ;
  std::vector<float> daqChannelWireGeomEndZ;
  std::vector<float> daqChannelWireGeomLength;
  std::vector<float> daqChannelSegmentStartZ;
  std::vector<float> daqChannelSegmentEndZ;
  std::vector<float> daqChannelSegmentStartY;
  std::vector<float> daqChannelSegmentEndY;
  std::vector<float> daqChannelSegmentLength;
  std::vector<float> daqChannelSegmentTension;
  std::vector<int>   daqChannelSegmentNumber;

  for (int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);

    if (segmentChannel+wireOffset == (int)thisHit->Channel() &&
        geom->View(thisHit->Channel()) == geom->View(thisHit->Channel()-wireOffset)){

      // information from the simulated geometry
      // wires on each plane have different start and end points:
      //
      // plane 0:
      //    -- start : 76.1 mm 
      //    -- end   : 6066.7 mm
      // plane 1:
      //    -- start : 6063.35 mm
      //    -- end   : 76.1 mm
      // plane 2:
      //    -- start : 6060 mm
      //    -- end   : 76.1 mm
      //
      // here to make things consistent later on, I'm flipping plane 0
      // and saying it starts at the end, and ends at the start
      //
      // now correct information from excel files by the correct amount

      double startPosition;
      if (thisHit->View() == 0)
        startPosition = 6066.7;
      if (thisHit->View() == 1)
        startPosition = 6063.35;
      if (thisHit->View() == 2)
        startPosition = 6060;
      
      segmentStartY = (startPosition)-segmentStartY;
      segmentEndY   = (startPosition)-segmentEndY;

      MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
        << "-- printing information from file "
        << "\n----starty : " << segmentStartY
        << "\n----startz : " << segmentStartZ
        << "\n----endy   : " << segmentEndY
        << "\n----endz   : " << segmentEndZ;

      std::vector<geo::WireID> wireIDs = geom->ChannelToWire(thisHit->Channel());

      MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
        << "-- found " << wireIDs.size() << " candidate wires corresponding to channel " << thisHit->Channel() << "\n"
        << "-- printing information from geometry";

      // define variables for the loop
      double distBetweenSegmentsY  = 99999;
      double closestWireGeomStartY = -1;
      double closestWireGeomEndY   = -1;
      double closestWireGeomStartZ = -1;
      double closestWireGeomLength = -1;
      double closestWireGeomEndZ   = -1;
      double segmentStartYSegment  = -1;
      double segmentStartZSegment  = -1;
      double segmentEndYSegment    = -1;
      double segmentEndZSegment    = -1;
      double segmentLength         = -1;

      // loop wire IDs and find the correct one

      for (size_t iwid = 0; iwid < wireIDs.size(); iwid++){

        geo::WireID thisWireID = wireIDs.at(iwid);

        geo::WireGeo const& thisWireGeo = geom->Wire(thisWireID);

        // convert to mm
        double thisWireGeoStartX = thisWireGeo.GetStart()[0]*10;
        double thisWireGeoStartY = thisWireGeo.GetStart()[1]*10;
        double thisWireGeoStartZ = thisWireGeo.GetStart()[2]*10;
        double thisWireGeoEndX   = thisWireGeo.GetEnd()[0]*10;
        double thisWireGeoEndY   = thisWireGeo.GetEnd()[1]*10;
        double thisWireGeoEndZ   = thisWireGeo.GetEnd()[2]*10;
        double thisWireGeoLength = thisWireGeo.Length();

        MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
          << "---- WireID "  << iwid
          << "\n"
          << "---- length: " << thisWireGeo.Length() << "\n"
          << "---- start: "  << thisWireGeoStartX
          << ", "            << thisWireGeoStartY
          << ", "            << thisWireGeoStartZ    << "\n"
          << "---- end: "    << thisWireGeoEndX
          << ", "            << thisWireGeoEndY
          << ", "            << thisWireGeoEndZ;

        double newDistY         = -1;
        float rotatedGeomStartY = -1;
        float rotatedGeomEndY   = -1;
        float rotatedGeomStartZ = -1;
        float rotatedGeomEndZ   = -1;

        // different rotations for U vs V
        if (thisHit->View() == 1){
          newDistY = std::abs(segmentStartY - thisWireGeoEndY) + std::abs(segmentEndY - thisWireGeoStartY);
          rotatedGeomEndY     = thisWireGeoStartY;
          rotatedGeomStartY   = thisWireGeoEndY;
          rotatedGeomEndZ     = thisWireGeoStartZ;
          rotatedGeomStartZ   = thisWireGeoEndZ;
        }
        else if (thisHit->View() == 0){
          newDistY = std::abs(segmentStartY - thisWireGeoStartY) + std::abs(segmentEndY - thisWireGeoEndY);
          rotatedGeomEndY     = thisWireGeoEndY;
          rotatedGeomStartY   = thisWireGeoStartY;
          rotatedGeomEndZ     = thisWireGeoEndZ;
          rotatedGeomStartZ   = thisWireGeoStartZ;
        }

        if (newDistY < distBetweenSegmentsY){
          distBetweenSegmentsY  = newDistY;
          closestWireGeomStartY = rotatedGeomStartY;
          closestWireGeomEndY   = rotatedGeomEndY;
          closestWireGeomStartZ = rotatedGeomStartZ;
          closestWireGeomEndZ   = rotatedGeomEndZ;
          closestWireGeomLength = thisWireGeoLength;
          segmentStartYSegment  = segmentStartY;
          segmentStartZSegment  = segmentStartZ;
          segmentEndYSegment    = segmentEndY;
          segmentEndZSegment    = segmentEndZ;
          segmentLength         = std::sqrt(std::pow(segmentEndZ-segmentStartZ,2) + std::pow(segmentEndY-segmentStartY,2));
        }

      }

      daqChannelWireGeomStartY . push_back(closestWireGeomStartY);
      daqChannelWireGeomEndY   . push_back(closestWireGeomEndY);
      daqChannelWireGeomStartZ . push_back(closestWireGeomStartZ);
      daqChannelWireGeomEndZ   . push_back(closestWireGeomEndZ);
      daqChannelWireGeomLength . push_back(closestWireGeomLength);
      daqChannelSegmentStartZ  . push_back(segmentStartZSegment);
      daqChannelSegmentEndZ    . push_back(segmentEndZSegment  );
      daqChannelSegmentStartY  . push_back(segmentStartYSegment);
      daqChannelSegmentEndY    . push_back(segmentEndYSegment  );
      daqChannelSegmentLength  . push_back(segmentLength);
      daqChannelSegmentNumber  . push_back(segmentNumber);
      daqChannelSegmentTension . push_back(segmentTension);
    }
  }

  MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
    << "-- found " << daqChannelSegmentTension.size() << " matched segments";

  if (daqChannelSegmentTension.size() == 0){

    MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
      << "-- that's a problem, dumping information: \n" 
      << "-- thisHit daq channel: " << thisHit->Channel() << "\n"
      << "-- thisHit view: " << thisHit->View();

    return thisTensionInformation; 

    //throw std::logic_error("Zero segments associated with DAQ channel.");

  }

  if ( thisHit->View() !=2){
    // get spacepoint z and y positions and APA wire x and y positions
    // remember that APA(x) is spacepoint(z)
    // nameing spacepoint_z to avoid confusion

    MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
      << "-- hit is on plane: " << thisHit->View()
      << "\n-- hit has channel: " << thisHit->Channel();

    // convert to mm
    float spacepoint_z      = thisSpacePoint->XYZ()[2]*10;
    float spacepoint_y      = thisSpacePoint->XYZ()[1]*10;
    float thisDist          = 9999;
    float thisTension       = 0;
    int   thisSegment       = 0;
    float thisGeomStartY    = 0;
    float thisGeomEndY      = 0;
    float thisGeomStartZ    = 0;
    float thisGeomEndZ      = 0;
    float thisGeomLength    = 0;
    float thisSegmentStartY = 0;
    float thisSegmentEndY   = 0;
    float thisSegmentStartZ = 0;
    float thisSegmentEndZ   = 0;
    float thisSegmentLength = 0;

    MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
      << "-- spacepoint_z: " << spacepoint_z
      << "\n-- spacepoint_y: " << spacepoint_y;

    // loop over wire segments and find the distance between the wire segment and
    // the spacepoint in order to find the segment closest to the hit

    for (size_t i = 0; i < daqChannelSegmentTension.size(); i++){

      MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
        << "-- iteration " << i;

      float line_y1          = daqChannelWireGeomStartY . at(i);
      float line_y2          = daqChannelWireGeomEndY   . at(i);
      float line_z1          = daqChannelWireGeomStartZ . at(i);
      float line_z2          = daqChannelWireGeomEndZ   . at(i);
      float daqTension       = daqChannelSegmentTension . at(i);
      float daqGeomStartY    = daqChannelWireGeomStartY . at(i);
      float daqGeomEndY      = daqChannelWireGeomEndY   . at(i);
      float daqGeomStartZ    = daqChannelWireGeomStartZ . at(i);
      float daqGeomEndZ      = daqChannelWireGeomEndZ   . at(i);
      float daqGeomLength    = daqChannelWireGeomLength . at(i);
      float daqSegmentStartY = daqChannelSegmentStartY  . at(i);
      float daqSegmentEndY   = daqChannelSegmentEndY    . at(i);
      float daqSegmentStartZ = daqChannelSegmentStartZ  . at(i);
      float daqSegmentEndZ   = daqChannelSegmentEndZ    . at(i);
      float daqSegmentLength = daqChannelSegmentLength  . at(i);
      int   daqSegment       = daqChannelSegmentNumber  . at(i);

      MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
        << "-- segment defined by (" 
        << line_z1 
        << ", " 
        << line_y1
        << ") (" 
        << line_z2 
        << ", " 
        << line_y2
        << ")";

      // distance between a point and a line defined
      // by two points
      // from wikipedia
      float calcDist = std::abs(((line_y2-line_y1)*spacepoint_z) - ((line_z2-line_z1)*spacepoint_y) + (line_z2*line_y1) - (line_y2*line_z1)) / std::sqrt(std::pow((line_y2-line_y1),2) + std::pow((line_z2-line_z1),2));

      MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
        << "-- calculated distance between spacepoint and segment: "
        << calcDist;

      if (calcDist < thisDist) {
        thisDist          = calcDist;
        thisTension       = daqTension;
        thisSegment       = daqSegment;
        thisGeomStartY    = daqGeomStartY;
        thisGeomEndY      = daqGeomEndY;
        thisGeomStartZ    = daqGeomStartZ;
        thisGeomEndZ      = daqGeomEndZ;
        thisGeomLength    = daqGeomLength;
        thisSegmentStartY = daqSegmentStartY;
        thisSegmentEndY   = daqSegmentEndY;
        thisSegmentStartZ = daqSegmentStartZ;
        thisSegmentEndZ   = daqSegmentEndZ;
        thisSegmentLength = daqSegmentLength;

        MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
          << "-- spacepoint is closer to new segment " << thisSegment;
      }
    }

    // now have 
    // the distance between the space point and the closest segment
    // the tension associated with that segment
    // the number associated with that segment

    thisTensionInformation.wireGeometryYStart = thisGeomStartY;
    thisTensionInformation.wireGeometryYEnd   = thisGeomEndY;
    thisTensionInformation.wireGeometryZStart = thisGeomStartZ;
    thisTensionInformation.wireGeometryZEnd   = thisGeomEndZ;
    thisTensionInformation.wireGeometryLength = thisGeomLength;
    thisTensionInformation.wireSegmentYStart  = thisSegmentStartY; 
    thisTensionInformation.wireSegmentYEnd    = thisSegmentEndY;
    thisTensionInformation.wireSegmentZStart  = thisSegmentStartZ;
    thisTensionInformation.wireSegmentZEnd    = thisSegmentEndZ;
    thisTensionInformation.wireSegmentLength  = thisSegmentLength;
    thisTensionInformation.distSPGeometry     = thisDist;
    thisTensionInformation.tension            = thisTension;
  }
  else {
    // in the case that we're on the X plane, then there's only one associated wire
    thisTensionInformation.distSPGeometry  = 0;
    thisTensionInformation.tension = daqChannelSegmentTension.at(0);
  }

  MF_LOG_DEBUG("pdsp::TensionAnalysis::getTensionInformation")
    << "-- calculated minimal distance: "   << thisTensionInformation.distSPGeometry
    << "\n-- tension:                     " << thisTensionInformation.tension;

  return thisTensionInformation;

}

pdsp::CalorimetryInformation pdsp::TensionAnalysis::getCalorimetryInformation(art::Ptr<recob::SpacePoint> thisSpacePoint, std::vector< art::Ptr<anab::Calorimetry> > calVec){

  pdsp::CalorimetryInformation thisCalorimetryInformation;

  double spx = thisSpacePoint->XYZ()[0];
  double spy = thisSpacePoint->XYZ()[1];
  double spz = thisSpacePoint->XYZ()[2];

  thisCalorimetryInformation.distanceFromSPToCalo = 99999;
  thisCalorimetryInformation.energyDeposition     = -1;
  thisCalorimetryInformation.chargeDeposition     = -1;
  thisCalorimetryInformation.residualRange        = -1;

  for (size_t i = 0; i < calVec.size(); i++){

    art::Ptr< anab::Calorimetry > thisCalorimetry = calVec.at(i);

    std::vector< anab::Point_t > calxyzVec = thisCalorimetry->XYZ();
    std::vector< float > dedxVec = thisCalorimetry->dEdx();
    std::vector< float > dqdxVec = thisCalorimetry->dQdx();
    std::vector< float > resRgVec = thisCalorimetry->ResidualRange();

    for (size_t j = 0; j < calxyzVec.size(); j++){

      double calx = calxyzVec.at(j).X();
      double caly = calxyzVec.at(j).Y();
      double calz = calxyzVec.at(j).Z();

      float newDist = std::sqrt(std::pow(spx-calx,2) + std::pow(spy-caly,2) + std::pow(spz-calz,2));

      if (newDist < thisCalorimetryInformation.distanceFromSPToCalo){
        thisCalorimetryInformation.distanceFromSPToCalo = newDist;
        thisCalorimetryInformation.energyDeposition     = dedxVec.at(j);
        thisCalorimetryInformation.chargeDeposition     = dqdxVec.at(j);
        thisCalorimetryInformation.residualRange        = resRgVec.at(j);        
      }
    }
  }

  return thisCalorimetryInformation;
}

DEFINE_ART_MODULE(pdsp::TensionAnalysis)
