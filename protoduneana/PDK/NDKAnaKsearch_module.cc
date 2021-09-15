//Module analyzer
//Ana module for nucleon decay and atmospheric analysis
//Ana TTree contains MC truth and reconstruction info
//ahiguera@central.uh.edu

//========================================================================
//                          Headers
//========================================================================
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TH1.h"

// LArSoft includes
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

//Headers added by Tyler
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larreco/ShowerFinder/ShowerReco3D/ShowerRecoAlg.h"
#include "lardata/Utilities/PxUtils.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larreco/ShowerFinder/ShowerReco3D/ShowerRecoAlgBase.h"
#include "larreco/ShowerFinder/ShowerReco3D/ShowerCalo.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/ArtDataHelper/MVAReader.h"
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#define MAX_TRACKS 40
#define MAX_GEN_PARTS 100
#define MAX_MC_PARTS 2500
#define MAX_MC_TRAJ_PTS 100
#define MAX_FLASHES 200
#define MAX_SHOWERS 100
#define MVA_LENGTH 4
#define MAX_HITS 1000
#define MAX_CALO_PTS 500

using std::cout;
using std::endl;
using std::map;
using std::max;
using std::set;
using std::size;
using std::sort;
using std::string;
using std::vector;
using fhicl::ParameterSet;
using art::fill_ptr_vector;
using art::FindManyP;
using art::Handle;
using art::Ptr;
using art::ServiceHandle;
using detinfo::DetectorClocksData;
using cheat::ParticleInventoryService;
using cheat::BackTrackerService;
using simb::MCNeutrino;
using simb::MCParticle;
using simb::MCTruth;
using sim::TrackIDE;
using recob::Hit;
using recob::Track;
using recob::Vertex;
using anab::Calorimetry;
using anab::ParticleID;

namespace DUNE {

  class NDKAna : public art::EDAnalyzer {
  public:

    explicit NDKAna(ParameterSet const& pset);
    virtual ~NDKAna();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void MCTruthInformation(DetectorClocksData const& clockData, const MCParticle *particle);
    void reconfigure(ParameterSet const& pset);
    void PIDAcal(vector<const Calorimetry*> cal, vector<double> &PIDA);
    void Process(const art::Event& evt, bool &isFiducial);
    void truthMatcher(DetectorClocksData const& clockData, 
      vector<Ptr<Hit>> all_hits, vector<Ptr<Hit>> track_hits,
      const MCParticle*& MCparticle, double &Efrac, double &Ecomplet);
    void truthMatcher2(DetectorClocksData const& clockData,
      vector<Ptr<Track>> &tracklist, FindManyP<Hit> &track_hits,
      map<int, int> &mc_id2idx );
    double truthLength( const MCParticle *MCparticle );
    bool insideFV(double vertex[4]);

   // art::InputTag fSimulationProducerLabel;
  private:
    void reset();

    //========================================================================
    //                      Configuration Parameters: Begin                       
    //========================================================================

    /*
    // the parameters we'll read from the .fcl
    string PandoraLabel = "pandora";
    art::InputTag fPandoraLabel{PandoraLabel};
      
    art::InputTag const fPandoraTrackModuleLabel = {
      fPandoraLabel.label() + "Track",
      fPandoraLabel.instance(),
      fPandoraLabel.process()};

    art::InputTag const fPandoraShowerModuleLabel = {
      fPandoraLabel.label() + "Shower",
      fPandoraLabel.instance(),
      fPandoraLabel.process()};

    art::InputTag const fPandoraCaloModuleLabel = {
      fPandoraLabel.label() + "calo",
      fPandoraLabel.instance(),
      fPandoraLabel.process()};

    art::InputTag const fPandoraPIDModuleLabel = {
      fPandoraLabel.label() + "pid",
      fPandoraLabel.instance(),
      fPandoraLabel.process()};

    */

    string fHitModuleLabel;
    string fTrackModuleLabel;
    string fPidModuleLabel;
    string fCaloModuleLabel;
    string fOpFlashModuleLabel;
    string fShowerModuleLabel;
    string fPointIDModuleLabel;
    string fNNetModuleLabel;
    string fPFParticleModuleLabel;
    string fMCgenieLabel;
    double fExponentConstant;
    double fPIDA_endPoint;
    double fMaxPIDAValue;
    double fMinPIDAValue;
    bool fSaveMCTree;

    // TH1F* h_MC_end_Momentum_Muons;
    // TH1F* h_MC_Momentum_Muons;    
    // TH1F* h_MC_end_Momentum;
    // TH1F* h_MC_Momentum;
    // TH1F* h_MC_end_Momentum_Kaons;
    // TH1F* h_MC_Momentum_Kaons;
    // TH1F* h_MC_length_Neutrinos;
    TH1F* h_MC_DaughterMuons;
    TH1F* h_MC_MotherKaons;
    TH1F* h_MC_PDG;
    TH1F* h_MC_PDGM;
    // TH1F* h_MC_length_Electrons;
    // TH1F* h_MC_length_Protons;
    TH1F* h_MC_length_Kaons;
    TH1F* h_MC_length_Muons;
    TH1F* h_MC_length_Pions;
    // TH1F* h_MC_length_NegPions;
    TH1F* h_MC_length;
    TH1F* h_track_length;
    TH1F* h_track_length_Walls_Kaons;
    TH1F* h_track_length_NonWalls_Kaons;
    TH1F* h_track_length_Walls_Muons;
    TH1F* h_track_length_NonWalls_Muons;
    TH1F* h_shortMuonsX;
    TH1F* h_shortMuonsY;    
    TH1F* h_shortMuonsZ;
    TH1F* h_track_length_Kaons;
    TH1F* h_track_length_Muons;
    TH1F* h_track_length_Pions;
    //========================================================================
    //                      Configuration Parameters: Emd
    //========================================================================
    
    TTree *fEventTree;

    // Event
    int Event;    //Event Number
    int Run;      //Run Number??
    int SubRun;   //SubRun Number??
    //========================================================================
    //                      Tree Branch Buffers: Begin
    //========================================================================



    //========================================================================
    //                           Truth Variables
    ////======================================================================

    //======= Neutrino True Information =========
     
    double MC_Ev;          //Neutrino Momentum
    int    MC_nuPDG;       //Neutrino flavor (check PDG for number)
    double MC_Q2;          //Momentum squared transfered in the neutrino interaction
    double MC_hit_nucleon; //The interaction was with a proton or neutrons?
    int    MC_cc;          //1 for CC neutrino interaction, 0 for NC interaction
      
    //======= Other Particle True Information =========
    int    MCgenie_npart;                            //Number of particles in each event
    int    MCgenie_id[MAX_GEN_PARTS];                //Track ID of particles
    int    MCgenie_pdg[MAX_GEN_PARTS];               //Particle ID of particle
    int    MCgenie_mother[MAX_GEN_PARTS];            //Is the particle a mother? yes = 0, No = Mother number
    int    MCgenie_statusCode[MAX_GEN_PARTS];        //Get's the status code retruned by GENIE/Geant4
    int    MCgenie_fate[MAX_GEN_PARTS];              //Rescatter?
    double MCgenie_startMomentum[MAX_GEN_PARTS][4];  //starting momentum
    double MCgenie_endMomentum[MAX_GEN_PARTS][4];    //ending momentum

    //======= Geant4 True Information =========
    double MC_vertex[4];                              //True Vertex Position
    // double MC_Vertex[4];
    int    MC_npart;                                  //Number of particles in each event
    int    MC_id[MAX_MC_PARTS];                       //Track ID of particles
    int    MC_pdg[MAX_MC_PARTS];                      //Particle ID of particles
    int    MC_mother[MAX_MC_PARTS];                   //Is the particle a mother?  yes = 0, No = Mother number
    double MC_startXYZT[MAX_MC_PARTS][4];             //Starting position of track
    double MC_endXYZT[MAX_MC_PARTS][4];               //Ending position of the track
    int    MC_npoints[MAX_MC_PARTS];                  //??
    double MC_xyz[MAX_MC_PARTS][MAX_MC_TRAJ_PTS][3];  //??
    double MC_startMomentum[MAX_MC_PARTS][4];         //Starting momentum of track
    double MC_endMomentum[MAX_MC_PARTS][4];           //Ending momentum of track
    double MC_truthlength[MAX_MC_PARTS];              //Length of track
    double MC_length_Kaons[MAX_MC_PARTS];
    double MC_length_Muons[MAX_MC_PARTS];
    double MC_length_Pions[MAX_MC_PARTS];
    double MC_length_KMP[MAX_MC_PARTS];
    // int MC_PDG[MAX_MC_PARTS];
    // int MC_PDGM[MAX_MC_PARTS];

    double MC_Prange[MAX_MC_PARTS];                   //Not Implemented??
    int    MC_statusCode[MAX_MC_PARTS];               //Get the status code teruned by GENIE/Geant4
    vector<string> MC_process;            //Process that created the particle, "Mother Particle" = "primary"
    vector<string> MC_Endprocess;         //Process that destroyed the particle

    // mc to reco associations

    //========================================================================
    //                     Pandora PFParticle Variables
    //========================================================================             
    int nPFParticle; //Number of PFParticles
    vector<int> PFP_PdgCode;  //A vector of PFParticle PDG codes
    vector<int> PFPTrkSize;
    vector<int> PFPNumber;
    vector<double> PFPTrkLengths;
    vector<double> PFPStartMomentum;
    vector<double> PFPEndMomentum;
    vector<TVector3> PFPStartPoint;
    vector<TVector3> PFPEndPoint;
    vector<TVector3> PFPStartDirection;
    vector<int> PFP_Parent;
    vector<int> PFParticleID;
    vector<int> nDaughters;
    vector<int> nTrk_Cal;
    //vector<bool> Track_isContained;
   
    int nPrimaries;

  //========================================================================
  //               Backtracking and Reconstructed Variables
  //======================================================================== 
    int    MC_nrecotracks[MAX_MC_PARTS];
    int    MC_reco_track_ids[MAX_MC_PARTS][MAX_TRACKS];
    int    MC_nhits_in_track[MAX_MC_PARTS][MAX_TRACKS];
    double MC_energy_in_track[MAX_MC_PARTS][MAX_TRACKS];


    int    n_vertices;                                              //Number of reconstructed vertices
    double vertex[MAX_TRACKS][4];                                   //Position of the Vertex??
    // double Vertex[Max_TRACKS][4];
    
    // int n_vertices_Kaons;
    double vertex_Kaons[MAX_TRACKS][4];
    // double Kaon_Tvertex[MAX_TRACKS][4]; 
    // double Kaon_Rvertex[MAX_TRACKS][4];  

    int    vtx_ID[MAX_TRACKS];                                      //Vertex ID
    int    n_decayVtx;                                              //Number of vertex from particle decays
    double decayVtx[MAX_TRACKS][3];                                 //Position of the decay vertex
    int    n_recoTracks;                                            //number of reconstructed tracks
    int    track_isContained[MAX_TRACKS];                           //number of reconstructed tracks fully contained within the detector
    int    track_ID[MAX_TRACKS];                                    //reconstructed track ID
    int    vtxID_trk[MAX_TRACKS][10];                               //Vertex Track ID
    double track_vtxDir[MAX_TRACKS][3];                             //Vertex direction of the reconstructed track
    double track_vtx[MAX_TRACKS][4];                                //Vertex position of the reconstructed track
    double track_end[MAX_TRACKS][4];                                //End of the reconstructed track
    double track_length[MAX_TRACKS];                                //Length of the reconstructed track
    
    double track_length_Kaons[MAX_TRACKS];   
    double track_length_Muons[MAX_TRACKS];
    double track_length_Pions[MAX_TRACKS];
    double track_length_KMP[MAX_TRACKS];

    double track_dir_vtx[MAX_TRACKS][4];                            //??
    double track_PIDA[MAX_TRACKS][3];                               //PIDA value of the track
    int    track_PID_pdg[MAX_TRACKS][3];                            //PDG from ParticleID algorithm
    double track_KE[MAX_TRACKS][3];                                 //Total Kinetic Energy
    int    track_bestplane[MAX_TRACKS];                             //An estimation of the best plane based on DOF????
    double track_Prange[MAX_TRACKS];                                //P = momentum range, Momentum calculated assuming the range from Track and a particle
    double track_Efrac[MAX_TRACKS];                                 //ratio between the amount of energy deposited by the particle that deposited more energy and the ttoal energy deposited by all particles
    double track_complet[MAX_TRACKS];                               //ratio between the amount of energy deposited by the particle that deposited more energy and the total energy deposited in the event
    int    track_mcID[MAX_TRACKS];                                  //True track ID used to compare with reoncstructed info 
    int    track_mcPDG[MAX_TRACKS];                                 //PDG particle used to compare with reconstructed info
    int    n_cal_points[MAX_TRACKS];                                //Number of points where depsoited energy in the best plane
    double track_dQ_dx[MAX_TRACKS][MAX_CALO_PTS];                   //Energy deposited with no correction
    double track_dE_dx[MAX_TRACKS][MAX_CALO_PTS];                   //Energy deposited
    double track_range[MAX_TRACKS][MAX_CALO_PTS];                   //Residual Range
    double track_pitch[MAX_TRACKS][MAX_CALO_PTS];                   //track pitch on the collection plane
    int    n_cal_points_byplane[MAX_TRACKS][3];                     //Number of points where deposited energy is in each plane
    double track_dQ_dx_byplane[MAX_TRACKS][3][MAX_CALO_PTS];        //Non-corrected energy deposited by plane
    double track_dE_dx_byplane[MAX_TRACKS][3][MAX_CALO_PTS];        //Energy deposited by plane
    double track_range_byplane[MAX_TRACKS][3][MAX_CALO_PTS];        //Residual range by plane
    double track_pitch_byplane[MAX_TRACKS][3][MAX_CALO_PTS];        //track pitch by plane
    double track_calo_xyz_byplane[MAX_TRACKS][3][MAX_CALO_PTS][3];  // storing xyz coordinates of all calo points

  //======= Reco Hit Information =========

    int     n_recoHits;
    int     hit_channel[MAX_HITS];
    int     hit_tpc[MAX_HITS];
    int     hit_plane[MAX_HITS];
    int     hit_wire[MAX_HITS];
    double  hit_peakT[MAX_HITS];
    double  hit_charge[MAX_HITS];
    double  hit_ph[MAX_HITS];
    double  hit_startT[MAX_HITS];
    double  hit_endT[MAX_HITS];
    double  hit_rms[MAX_HITS];
    double  hit_electrons[MAX_HITS];

    double Em_ch;
    double Em_e;
    double trk_e;
    double Emichel_e;

  //======= Reco Shower Information =========
    int    n_recoShowers;
    double sh_direction_X[MAX_SHOWERS];
    double sh_direction_Y[MAX_SHOWERS];
    double sh_direction_Z[MAX_SHOWERS];
    double sh_start_X[MAX_SHOWERS];
    double sh_start_Y[MAX_SHOWERS];
    double sh_start_Z[MAX_SHOWERS];
    double sh_energy[MAX_SHOWERS][3];
    double sh_MIPenergy[MAX_SHOWERS][3];
    double sh_dEdx[MAX_SHOWERS][3];
    int    sh_bestplane[MAX_SHOWERS];
    double sh_length[MAX_SHOWERS];

  //======= Reco Light Information? =========

    int    n_flashes;
    double flash_time[MAX_FLASHES];
    double flash_pe[MAX_FLASHES];
    double flash_ycenter[MAX_FLASHES];
    double flash_zcenter[MAX_FLASHES];
    double flash_ywidth[MAX_FLASHES];
    double flash_zwidth[MAX_FLASHES];
    double flash_timewidth[MAX_FLASHES];
    double flash_abstime[MAX_FLASHES];
    int    flash_frame[MAX_FLASHES];
    double flash_PE_ndk[MAX_FLASHES];
    double flash_PE_Ar39[MAX_FLASHES];

  //======= Far Detector Geometry =========

    double fFidVolCutX;
    double fFidVolCutY;
    double fFidVolCutZ;

    double fFidVolXmin;
    double fFidVolXmax;
    double fFidVolYmin;
    double fFidVolYmax;
    double fFidVolZmin;
    double fFidVolZmax;

    double fPidValue;
    unsigned int fView;

    bool fSaveHits;

    //========================================================================
    //                      Tree Branch Buffers: End
    //========================================================================


    //detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    DetectorClocksData const ts = ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    detinfo::DetectorPropertiesData const detprop = ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(ts);
    double XDriftVelocity = detprop.DriftVelocity()*1e-3; //cm/ns
    double WindowSize     = detprop.NumberTimeSamples() * ts.TPCClock().TickPeriod() * 1e3;
    
    ServiceHandle<geo::Geometry> geom;
    calo::CalorimetryAlg fCalorimetryAlg;

    protoana::ProtoDUNETruthUtils fTruthUtil; // utility to help with connecting MC truth and reco objects

  }; // class NDKAna

  //=======================================================================
  //PFPParticle Data Structure?  Added by Tyler D. Stokes
  //=======================================================================
  /* 
  template <typename T>
  using PFParticleData_t = vector<T>;
  template <typename T>
  using DaughterData_t = vector<BoxedArray<T<kMaxNDaughtersPerPFP]>>;
  template <typename T>
  using ClusterData_t = vector<BoxedArray<T[kMaxNClustersPerPFP]>>;

  size_t MaxPFParticles; // maximum number of storable PFParticles

  // @{
  // @name Branch data structures

   Short_t                   nPFParticles;     ///< the total number of PFParticles
   PFParticleData_t<Short_t> pfp_selfID;       ///< the PFParticles' own IDs
   PFParticleData_t<Short_t> pfp_isPrimary;    ///< whether the PFParticle is a primary particle
   PFParticleData_t<Short_t> pfp_numDaughters; ///< the number of daughters belonging to this PFParticle
   DaughterData_t<Short_t>   pfp_daughterIDs;  ///< the IDs of the daughter PFParticles
   PFParticleData_t<Short_t> pfp_parentID;     ///< the ID of this PFParticle's immediate parent
   PFParticleData_t<Short_t> pfp_vertexID;     ///< the ID of the vertex belonging to this PFParticle
   PFParticleData_t<Short_t> pfp_isShower;     ///< whether this PFParticle corresponds to a shower
   PFParticleData_t<Short_t> pfp_isTrack;      ///< whether this PFParticle corresponds to a track
   PFParticleData_t<Short_t> pfp_trackID;      ///< the ID of the track object corresponding to this PFParticle, if !isShower
   PFParticleData_t<Short_t> pfp_showerID;     ///< the ID of the shower object corresponding to this PFParticle, if isShower
   PFParticleData_t<Short_t> pfp_isNeutrino;   ///< whether this PFParticle is a neutrino
   PFParticleData_t<Int_t>   pfp_pdgCode;      ///< the preliminary estimate of the PFParticle type using the PDG code
   PFParticleData_t<Short_t> pfp_numClusters;  ///< the number of associated clusters
   ClusterData_t<Short_t>    pfp_clusterIDs;   ///< the IDs of any associated clusters
   Short_t                   pfp_numNeutrinos; ///< the number of reconstructed neutrinos
   Short_t pfp_neutrinoIDs[kMaxNPFPNeutrinos]; ///< the PFParticle IDs of the neutrinos
  //@}

  //Creates a PFParticle data structure allowing up to maxPFParticles PFParticles
  PFParticleDataStruct(size_t maxPFParticles = 0):
    MaxPFParticles(maxPFParticles) { Clear(); }

  void Clear();
  void SetMaxPFParticles(size_t maxPFParticles)
    { MaxPFParticles = maxPFParticles; Resize(MaxPFParticles); }
  void Resize(size_t numPFParticles);
  void SetAddresses(TTree* pTree);

  size_t GetMaxPPFParticles() const { return MaxPFParticles; }
  size_t GetMaxDaughtersPerPFParticle(int //iPFParticle // = 0) const
    { return (size_t_ kMaxNDaughterPerPFP; }
  size_t GetMaxClusterPerPFParticle(int //iPFParticle // = 0) const
    { return (size_t) kMaxNClustersPerPFP; }

  }; //class PFParticleDataStruct                                                                                                      
   */


  NDKAna::NDKAna(ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet), fCalorimetryAlg(parameterSet.get< ParameterSet >("CalorimetryAlg"))
  {
  // fSimulationProducerLabel=parameterSet.get<string>("SimulationLabel");

    reconfigure(parameterSet);

  }

  //========================================================================
  NDKAna::~NDKAna() {
    //destructor
  }
  //========================================================================
  void NDKAna::reconfigure(ParameterSet const& p)
  {
    fHitModuleLabel      = p.get<string>("HitModuleLabel");
    fTrackModuleLabel    = p.get<string>("TrackModuleLabel");
    fPidModuleLabel      = p.get<string>("PidModuleLabel");
    fCaloModuleLabel     = p.get<string>("CaloModuleLabel");
    fOpFlashModuleLabel  = p.get<string>("OpFlashModuleLabel");
    fShowerModuleLabel   = p.get<string>("ShowerModuleLabel");
    fPointIDModuleLabel  = p.get<string>("PointIDModuleLabel");
    fNNetModuleLabel     = p.get<string>("NNetModuleLabel");
    fMCgenieLabel        = p.get<string>("MCgenieLabel");
    fSaveHits            = p.get<bool>("SaveHits");
    fExponentConstant    = p.get<double>("ExponentConstant");
    fMaxPIDAValue        = p.get<double>("MaxPIDAValue");
    fMinPIDAValue        = p.get<double>("MinPIDAValue");
    fPIDA_endPoint       = p.get<double>("PIDA_endPoint");
    fView                = p.get<double>("View");
    fPidValue            = p.get<double>("PidValue");
    fFidVolCutX          = p.get<double>("FidVolCutX");
    fFidVolCutY          = p.get<double>("FidVolCutY");
    fFidVolCutZ          = p.get<double>("FidVolCutZ");
  }
  //========================================================================
  void NDKAna::beginJob() {

    cout<<"job begin..."<<endl;
    // Get geometry.
    auto const* geo = lar::providerFrom<geo::Geometry>();
    // Define histogram boundaries (cm).
    // For now only draw cryostat=0.
    double minx = 1e9;
    double maxx = -1e9;
    double miny = 1e9;
    double maxy = -1e9;
    double minz = 1e9;
    double maxz = -1e9;
    for (size_t i = 0; i<geo->NTPC(); ++i){
      double local[3] = {0.,0.,0.};
      double world[3] = {0.,0.,0.};
      const geo::TPCGeo &tpc = geo->TPC(i);
      tpc.LocalToWorld(local,world);
      if (minx>world[0]-geo->DetHalfWidth(i))
        minx = world[0]-geo->DetHalfWidth(i);
      if (maxx<world[0]+geo->DetHalfWidth(i))
        maxx = world[0]+geo->DetHalfWidth(i);
      if (miny>world[1]-geo->DetHalfHeight(i))
        miny = world[1]-geo->DetHalfHeight(i);
      if (maxy<world[1]+geo->DetHalfHeight(i))
        maxy = world[1]+geo->DetHalfHeight(i);
      if (minz>world[2]-geo->DetLength(i)/2.)
        minz = world[2]-geo->DetLength(i)/2.;
      if (maxz<world[2]+geo->DetLength(i)/2.)
        maxz = world[2]+geo->DetLength(i)/2.;
    }

    fFidVolXmin = minx + fFidVolCutX;
    fFidVolXmax = maxx - fFidVolCutX;
    fFidVolYmin = miny + fFidVolCutY;
    fFidVolYmax = maxy - fFidVolCutY;
    fFidVolZmin = minz + fFidVolCutZ;
    fFidVolZmax = maxz - fFidVolCutZ;

    cout<<"Fiducial volume:"<<"\n"
       <<fFidVolXmin<<"\t< x <\t"<<fFidVolXmax<<"\n"
       <<fFidVolYmin<<"\t< y <\t"<<fFidVolYmax<<"\n"
       <<fFidVolZmin<<"\t< z <\t"<<fFidVolZmax<<"\n";

    ServiceHandle<art::TFileService> tfs;

    h_MC_MotherKaons             = tfs->make<TH1F>("h_MC_MotherKaons" , "Muon Mother Kaon Track Length", 80, 1, 100);
    h_MC_DaughterMuons           = tfs->make<TH1F>("h_MC_DaughterKaons" , "Kaon Muon Daughters Track Length" , 80, 1, 100);
    //h_MC_end_Momentum          = tfs->make<TH1F>("h_MC_end_Momentum", "True End Momentum", 100, 0, 500);
    //h_MC_Momentum        = tfs->make<TH1F>("h_MC_Momentum", "True Momentum", 100, 0, 500);
    //h_MC_end_Momentum_Muons    = tfs->make<TH1F>("h_MC_end_Momentum", "True End Momentum", 100, 0, 500);
    //h_MC_Momentum_Muons  = tfs->make<TH1F>("h_MC_Momentum", "True Muon Momentum", 100, 0, 500);
    //h_MC_end_Momentum_Kaons    = tfs->make<TH1F>("h_MC_end_Momentum", "True End Momentum", 100, 0, 500);
    //h_MC_Momentum_Kaons  = tfs->make<TH1F>("h_MC_Momentum", "True Kaon Momentum", 100, 0, 500);
    h_MC_PDGM                  = tfs->make<TH1F>("h_MC_PDGM","True Numbers Non-Decays", 100, -20, 350);
    h_MC_PDG                   = tfs->make<TH1F>("h_MC_PDG", "True Decays Numbers", 100, -20, 350);  
    h_MC_length                = tfs->make<TH1F>("h_MC_length", "True Track Length", 80, 1, 100);
    //h_MC_length_Neutrinos      = tfs->make<TH1F>("h_MC_length_Neutrinos", "True Track Length Neutrinos", 80, 1, 100);
    //h_MC_length_Electrons      = tfs->make<TH1F>("h_MC_length_Electrons" , "True Track Length Electrons", 80, 1, 100);
    //h_MC_length_Protons        = tfs->make<TH1F>("h_MC_length_Protons" , "True Track Length Protons", 80, 1, 100);
    h_MC_length_Kaons          = tfs->make<TH1F>("h_MC_length_Kaons" , "True Track Length Kaons", 80, 1, 100);
    h_MC_length_Muons          = tfs->make<TH1F>("h_MC_length_Muons" , "True Track Length Muons", 80, 1, 100);
    h_MC_length_Pions          = tfs->make<TH1F>("h_MC_length_Pions" , "True Track Length Pions", 80, 1, 100);
    //h_MC_length_NegPions       = tfs->make<TH1F>("h_MC_length_NegPions" , "Truth Track Length Neg. Pins", 80, 1, 100);
    h_track_length             = tfs->make<TH1F>("h_track_length", "Reco Track Length", 80, 1, 80);
    h_track_length_Kaons       = tfs->make<TH1F>("h_track_length_Kaons", "Reco Track Length Kaons", 80, 1, 100);
    h_track_length_Muons       = tfs->make<TH1F>("h_track_length_Muons", "Reco Track Length Muons", 80, 1, 100);
    h_track_length_Pions       = tfs->make<TH1F>("h_track_length_Pions", "Reco Track Length Pions", 80, 1, 100);
    h_track_length_Walls_Muons = tfs->make<TH1F>("h_track_length_Walls_Muons" ,"Reco Track Length Edge Muons", 80, 1, 100);
    h_track_length_Walls_Kaons = tfs->make<TH1F>("h_track_length_Walls_Kaons" ,"Reco Track Length Edge Kaons", 80, 1, 100);
    h_track_length_NonWalls_Muons = tfs->make<TH1F>("h_track_length_NonWalls_Muons" ,"Reco Track LengthGood Muons", 80, 1, 100);
    h_track_length_NonWalls_Kaons = tfs->make<TH1F>("h_track_length_NonWalls_Kaons" ,"Reco Track Length Good Kaons", 80, 1, 100);
    h_shortMuonsX = tfs->make<TH1F>("h_shortMuonsX", "Reco Short Muons X" , 80, -360, 360);
    h_shortMuonsY = tfs->make<TH1F>("h_shortMuonsY", "Reco Short Muons Y" , 80, -600, 600);
    h_shortMuonsZ = tfs->make<TH1F>("h_shortMuonsZ", "Reco Short Muons Z" , 80, 0 , 1400);

    fEventTree = tfs->make<TTree>("Event", "Event Tree from Sim & Reco");

    fEventTree->Branch("eventNo", &Event);
    fEventTree->Branch("runNo", &Run);
    fEventTree->Branch("subRunNo", &SubRun);

    //GENIE info & list of particles
    fEventTree->Branch("MC_Ev", &MC_Ev);
    fEventTree->Branch("MC_cc", &MC_cc);
    fEventTree->Branch("MC_Q2", &MC_Q2);
    fEventTree->Branch("MC_nuPDG", &MC_nuPDG);
    fEventTree->Branch("MC_hit_nucleon", &MC_hit_nucleon);
    fEventTree->Branch("mcgenie_npart", &MCgenie_npart);  // number of particles
    fEventTree->Branch("mcgenie_id", &MCgenie_id, "mcgenie_id[mcgenie_npart]/I");
    fEventTree->Branch("mcgenie_fate", &MCgenie_fate, "mcgenie_fate[mcgenie_npart]/I");
    fEventTree->Branch("mcgenie_statusCode", &MCgenie_statusCode, "mcgenie_statusCode[mcgenie_npart]/I");
    fEventTree->Branch("mcgenie_pdg", &MCgenie_pdg, "mcgenie_pdg[mcgenie_npart]/I");
    fEventTree->Branch("mcgenie_mother", &MCgenie_mother, "mcgenie_mother[mcgenie_npart]/I");
    fEventTree->Branch("mcgenie_startMomentum", &MCgenie_startMomentum, "mcgenie_startMomentum[mcgenie_npart][4]/D");
    fEventTree->Branch("mcgenie_endMomentum", &MCgenie_endMomentum, "mcgenie_endMomentum[mcgenie_npart][4]/D");
   
    // Geant4 list of particles
    // fEventTree->Branch("mc_vertex", MC_vertex, "mc_vertex[4]/D");
    
    fEventTree->Branch("mc_vertex", MC_vertex, "mc_vertex[4]/D");

    fEventTree->Branch("mc_npart", &MC_npart);  // number of particles
    fEventTree->Branch("mc_id", MC_id, "mc_id[mc_npart]/I");
    fEventTree->Branch("mc_pdg", MC_pdg, "mc_pdg[mc_npart]/I");
    fEventTree->Branch("mc_statusCode", MC_statusCode, "mc_statusCode[mc_npart]/I");
    fEventTree->Branch("mc_mother", MC_mother, "mc_mother[mc_npart]/I");
    fEventTree->Branch("mc_startXYZT", MC_startXYZT, "mc_startXYZT[mc_npart][4]/D");
    fEventTree->Branch("mc_endXYZT", MC_endXYZT, "mc_endXYZT[mc_npart][4]/D");
    fEventTree->Branch("mc_npoints", MC_npoints, "mc_npoints[mc_npart]/I");
    fEventTree->Branch("mc_xyz", MC_xyz, "mc_xyz[mc_npart][100][3]/D");
    fEventTree->Branch("mc_startMomentum", MC_startMomentum, "mc_startMomentum[mc_npart][4]/D");
    fEventTree->Branch("mc_endMomentum", MC_endMomentum, "mc_endMomentum[mc_npart][4]/D");
    fEventTree->Branch("mc_Prange", MC_Prange, "mc_Prange[mc_npart]/D");
    fEventTree->Branch("mc_truthlength", MC_truthlength, "mc_truthlength[mc_npart]/D");
    fEventTree->Branch("MC_length_Kaons", MC_length_Kaons, "MC_length_Kaons[mc_npart]/D");
    fEventTree->Branch("MC_length_Muons", MC_length_Muons, "MC_length_Muons[mc_npart]/D");
    fEventTree->Branch("MC_length_Pions", MC_length_Pions, "MC_length_Pions[mc_npart]/D");
    fEventTree->Branch("MC_length_KMP", MC_length_KMP,     "MC_length_KMP[mc_npart]/D");

    fEventTree->Branch("mc_process", &MC_process);
    fEventTree->Branch("mc_Endprocess", &MC_Endprocess);
    fEventTree->Branch("mc_nrecotracks", MC_nrecotracks, "mc_nrecotracks[mc_npart]/I");
    fEventTree->Branch("mc_reco_track_ids", MC_reco_track_ids, "mc_reco_track_ids[mc_npart][40]/I");
    fEventTree->Branch("mc_nhits_in_track", MC_nhits_in_track, "mc_nhits_in_track[mc_npart][40]/I");
    fEventTree->Branch("mc_energy_in_track", MC_energy_in_track, "mc_energy_in_track[mc_npart][40]/D");
    fEventTree->Branch("n_vertices", &n_vertices);
    fEventTree->Branch("vertex", vertex,"vertex[n_vertices][4]/D");
    fEventTree->Branch("vtx_ID", vtx_ID,"vtx_ID[n_vertices]/I");
    fEventTree->Branch("n_reco_tracks", &n_recoTracks);
    fEventTree->Branch("n_decayVtx", &n_decayVtx);
    fEventTree->Branch("decayVtx", decayVtx,"decayVtx[n_decayVtx][3]/D");  //vertices found using decayID point Alg
    fEventTree->Branch("vtxID_trk", vtxID_trk,"vtxID_trk[n_reco_tracks][10]/I"); //track-vertex association
    fEventTree->Branch("track_vtx", track_vtx,"track_vtx[n_reco_tracks][4]/D");
    fEventTree->Branch("track_vtxDir", track_vtxDir,"track_vtxDir[n_reco_tracks][3]/D");
    fEventTree->Branch("track_end", track_end,"track_end[n_reco_tracks][4]/D");
    fEventTree->Branch("track_isContained", track_isContained,"track_isContained[n_reco_tracks]/I");
    fEventTree->Branch("track_ID", track_ID,"track_ID[n_reco_tracks]/I");
    
    fEventTree->Branch("track_length_Kaons", track_length_Kaons,"track_length_Kaons[n_reco_tracks][3]/D");
    fEventTree->Branch("track_length_Pions", track_length_Pions,"track_length_Pions[n_reco_tracks][3]/D");
    fEventTree->Branch("track_length_Muons", track_length_Muons,"track_length_Muons[n_reco_tracks][3]/D");
    fEventTree->Branch("track_length_KMP", track_length_KMP, "track_length_KMP[n_reco_tracks][3]/D"); 
   
    fEventTree->Branch("track_length", track_length, "track_length[n_reco_tracks][3]/D");
    fEventTree->Branch("track_PIDA", track_PIDA,"track_PIDA[n_reco_tracks][3]/D");
    // fEventTree->Branch("track_PID_pdg", track_PID_pdg,"track_PID_pdg[n_reco_tracks][3]/I");
    fEventTree->Branch("track_PID_pdg", track_PID_pdg, "trck_PID_pdg[n_reco_tracks][3]/I");
    fEventTree->Branch("track_KE", track_KE,"track_KE[n_reco_tracks][3]/D");
    fEventTree->Branch("track_Prange", track_Prange,"track_Prange[n_reco_tracks]/D");
    fEventTree->Branch("track_bestplane", track_bestplane,"track_bestplane[n_reco_tracks]/I");
    fEventTree->Branch("n_cal_points", n_cal_points,"n_cal_points[n_reco_tracks]/I");
    fEventTree->Branch("track_dQ_dx", track_dQ_dx,"track_dQ_dx[n_reco_tracks][500]/D");
    fEventTree->Branch("track_dE_dx", track_dE_dx,"track_dE_dx[n_reco_tracks][500]/D");
    fEventTree->Branch("track_range", track_range,"track_range[n_reco_tracks][500]/D");
    fEventTree->Branch("track_pitch", track_pitch,"track_pitch[n_reco_tracks][500]/D");

    fEventTree->Branch("n_cal_points_byplane", n_cal_points_byplane,"n_cal_points_byplane[n_reco_tracks][3]/I");
    fEventTree->Branch("track_dQ_dx_byplane", track_dQ_dx_byplane,"track_dQ_dx_byplane[n_reco_tracks][500][3]/D");
    fEventTree->Branch("track_dE_dx_byplane", track_dE_dx_byplane,"track_dE_dx_byplane[n_reco_tracks][500][3]/D");
    fEventTree->Branch("track_range_byplane", track_range_byplane,"track_range_byplane[n_reco_tracks][500][3]/D");
    fEventTree->Branch("track_pitch_byplane", track_pitch_byplane,"track_pitch_byplane[n_reco_tracks][500][3]/D");
    fEventTree->Branch("track_calo_xyz_byplane", track_calo_xyz_byplane,"track_calo_xyz_byplane[n_reco_tracks][3][500][3]/D");

    fEventTree->Branch("track_complet", track_complet,"track_complet[n_reco_tracks]/D");  //track quality variable (completeness)
    fEventTree->Branch("track_Efrac", track_Efrac,"track_Efrac[n_reco_tracks]/D");        //track quality variable (purity)
    fEventTree->Branch("track_mcID", track_mcID,"track_mcID[n_reco_tracks]/I");           //true MC ID for a given track
    fEventTree->Branch("track_mcPDG", track_mcPDG,"track_mcPDG[n_reco_tracks]/I");        //true MC PDG for a given track

    fEventTree->Branch("Em_ch", &Em_ch);
    fEventTree->Branch("Em_e", &Em_e);
    fEventTree->Branch("trk_e", &trk_e);
    fEventTree->Branch("Emichel_e", &Emichel_e);
    if (fSaveHits){
      fEventTree->Branch("n_recoHits", &n_recoHits);
      fEventTree->Branch("hit_channel", &hit_channel, "hit_channel[n_recoHits]/I");
      fEventTree->Branch("hit_tpc", &hit_tpc, "hit_tpc[n_recoHits]/I");
      fEventTree->Branch("hit_plane", &hit_plane, "hit_plane[n_recoHits]/I");
      fEventTree->Branch("hit_wire", &hit_wire, "hit_wire[n_recoHits]/I");
      fEventTree->Branch("hit_peakT", &hit_peakT, "hit_peakT[n_recoHits]/D");
      fEventTree->Branch("hit_charge", &hit_charge, "hit_charge[n_recoHits]/D");
      fEventTree->Branch("hit_ph", &hit_ph, "hit_ph[n_recoHits]/D");
      fEventTree->Branch("hit_endT", &hit_endT, "hit_endT[n_recoHits]/D");
      fEventTree->Branch("hit_rms", &hit_rms, "hit_rms[n_recoHits]/D");
      fEventTree->Branch("hit_electrons", &hit_electrons, "hit_electrons[n_recoHits]/D");
    }

    fEventTree->Branch("n_showers", &n_recoShowers);
    fEventTree->Branch("sh_direction_X", &sh_direction_X, "sh_direction_X[n_showers]/D");
    fEventTree->Branch("sh_direction_Y", &sh_direction_Y, "sh_direction_Y[n_showers]/D");
    fEventTree->Branch("sh_direction_Z", &sh_direction_Z, "sh_direction_Z[n_showers]/D");
    fEventTree->Branch("sh_start_X", &sh_start_X, "sh_start_X[n_showers]/D");
    fEventTree->Branch("sh_start_Y", &sh_start_Y, "sh_start_Y[n_showers]/D");
    fEventTree->Branch("sh_start_Z", &sh_start_Z, "sh_start_Z[n_showers]/D");
    fEventTree->Branch("sh_energy", &sh_energy, "sh_energy[n_showers][3]/D");
    fEventTree->Branch("sh_MIPenergy", &sh_MIPenergy, "sh_MIPenergy[n_showers][3]/D");
    fEventTree->Branch("sh_dEdx", &sh_dEdx, "sh_dEdx[n_showers][3]/D");
    fEventTree->Branch("sh_bestplane", &sh_bestplane, "sh_bestplane[n_showers]/I");
    fEventTree->Branch("sh_length", &sh_length, "sh_length[n_showers]/D");

    fEventTree->Branch("n_flashes", &n_flashes);
    fEventTree->Branch("flash_time", &flash_time,"flash_time[n_flashes]/D");
    fEventTree->Branch("flash_pe", &flash_pe,"flash_pe[n_flashes]/D");
    fEventTree->Branch("flash_ycenter", &flash_ycenter,"flash_ycenter[n_flashes]/D");
    fEventTree->Branch("flash_zcenter", &flash_zcenter,"flash_zcenter[n_flashes]/D");
    fEventTree->Branch("flash_ywidth", &flash_ywidth,"flash_ywidth[n_flashes]/D");
    fEventTree->Branch("flash_zwidth", &flash_zwidth,"flash_zwidth[n_flashes]/D");
    fEventTree->Branch("flash_timewidth", &flash_timewidth, "flash_timewidth[n_flashes]/D");
    fEventTree->Branch("flash_abstime", &flash_abstime, "flash_abstime[n_flashes]/D");
    fEventTree->Branch("flash_frame", &flash_frame, "flash_frame[n_flashes]/I");

    fEventTree->Branch("flash_PE_ndk", &flash_PE_ndk, "flash_PE_ndk[n_flashes]/D");
    fEventTree->Branch("flash_PE_Ar39", &flash_PE_Ar39, "flash_PE_Ar39[n_flashes]/D");

  }
  //========================================================================
  void NDKAna::endJob() {}
  //========================================================================
  void NDKAna::beginRun(const art::Run& /*run*/){
    mf::LogInfo("NDKAna") << "begin run..." << endl;
  }
  //========================================================================
  void NDKAna::analyze( const art::Event& event ){
    if (event.isRealData()) return;
    cout << "NDKAna::analyze()" << endl;

    reset(); // reset data to be stored in the tree

    // cout<<"NDKAna::analyze(): after reset"<<endl;

    Event  = event.id().event();
    Run    = event.run();
    SubRun = event.subRun();
    bool isFiducial = false;
    Process(event, isFiducial);

    // cout<<"NDKAna::analyze(): after processing TT"<<endl;

    if(isFiducial) fEventTree->Fill();
    // fEventTree->Fill();
    // cout<<"NDKAna::analyze(): after filling"<<endl;
  } // function NDAna::analyze

  //========================================================================
  void NDKAna::Process( const art::Event& event, bool &isFiducial){
    //int idbg = 0;

    //save GENIE stuf for atmos
    Handle<vector<MCTruth>> MCtruthHandle;
    vector<Ptr<MCTruth>> MCtruthlist;
    if(event.getByLabel(fMCgenieLabel, MCtruthHandle))
      fill_ptr_vector(MCtruthlist, MCtruthHandle);

    // cout<<"dbg a: "<<idbg++<<endl;

    // For now assume that there is only one neutrino interaction...
    Ptr<MCTruth> MCtruth;
    // cout<<"MCtruthlist Size" << MCtruthlist.size() << endl;
    if (MCtruthlist.size() > 0) {
      MCtruth = MCtruthlist[0];
      if (MCtruth->NeutrinoSet()) {
        MCNeutrino nu = MCtruth->GetNeutrino();
        if( nu.CCNC() == 0 ) MC_cc = 1;
        else if ( nu.CCNC() == 1 ) MC_cc = 0;
        MCParticle neutrino = nu.Nu();
        MC_nuPDG = nu.Nu().PdgCode();
        const TLorentzVector& nu_momentum = nu.Nu().Momentum(0);
        double MC_incoming_P[4];
        nu_momentum.GetXYZT(MC_incoming_P);
        MC_Ev = nu_momentum[3];
        MC_Q2 = nu.QSqr();
        MC_hit_nucleon = nu.HitNuc();
      }
      MCgenie_npart = MCtruth->NParticles();
      if (MCgenie_npart > MAX_GEN_PARTS) MCgenie_npart = MAX_GEN_PARTS;
      for (int i =0; i<MCgenie_npart; ++i) {
         MCParticle particle = MCtruth->GetParticle(i);
         MCgenie_id[i] = particle.TrackId();
         MCgenie_pdg[i] = particle.PdgCode();
         MCgenie_mother[i] = particle.Mother();
         MCgenie_statusCode[i] =particle.StatusCode();
         MCgenie_fate[i] = particle.Rescatter();
         const TLorentzVector& momentumStart = particle.Momentum(0);
         const TLorentzVector& momentumEnd   = particle.EndMomentum();
         //!Save the true vertex as the vertex using primaries ...hmmm do you have another suggestion?
         momentumStart.GetXYZT(MCgenie_startMomentum[i]);
         momentumEnd.GetXYZT(MCgenie_endMomentum[i]);
      }
    } // if MCTruth vector not empty

    // ServiceHandle<BackTrackerService> bt;
    ServiceHandle<ParticleInventoryService> part_inv;
    const sim::ParticleList& plist = part_inv->ParticleList();

    auto particleHandle = event.getValidHandle<vector<MCParticle>> ("largeant");
    
    map< int, const MCParticle* > particleMap;
    int fSimTrackID =0.;
    for(auto const& particle : (*particleHandle))
    {
      fSimTrackID = particle.TrackId();
      particleMap[fSimTrackID] = &particle;
    }

    map<int, int> mc_id2idx;

    //cout<<"dbg ab: "<<idbg++<<endl;

    int idxPart = 0; // particle index
    MC_npart = plist.size();
    if (MC_npart > MAX_MC_PARTS) return;
    for (auto const& [i, particle] : plist) {
    // for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar) {
      // particle = ipar->second;

      mc_id2idx[particle->TrackId()] = idxPart;

      MC_id[idxPart] = particle->TrackId();
      // MC_pdg[i] = particle->PdgCode();

      // if (particle->PdgCode() < 2500) {
      //   MC_pdg[i] = particle->PdgCode();
      // } else MC_pdg[i] = 2450;

      // if (particle->PdgCode() == 13 || particle->PdgCode() == 12 || particle->PdgCode() == 211 || particle->PdgCode() == 321 || particle->PdgCode() == 2212) {
      // if (particle->PdgCode() == 211 ||particle->PdgCode() == 11 || particle->PdgCode() == 12 || particle->PdgCode() == -14 || particle->PdgCode() ==-13 || particle->PdgCode() == 14 || particle->PdgCode() == 321) {
      MC_pdg[idxPart] = particle->PdgCode();
      //   MC_Momentum[i] = particle->Momentum().E();
      // } else MC_pdg[i] = -23;

      MC_mother[idxPart] = particle->Mother();
      MC_process.push_back(particle->Process());
      MC_Endprocess.push_back(particle->EndProcess());
      MC_statusCode[idxPart] = particle->StatusCode();
      const TLorentzVector& positionStart = particle->Position(0);
      const TLorentzVector& positionEnd = particle->EndPosition();
      const TLorentzVector& momentumStart = particle->Momentum(0);
      const TLorentzVector& momentumEnd = particle->EndMomentum();
      // Save the true vertex as the vertex using primaries ...hmmm do you have another suggestion?
      if( particle->Mother() == 0 ) positionStart.GetXYZT(MC_vertex);
      positionStart.GetXYZT(MC_startXYZT[idxPart]);
      positionEnd.GetXYZT(MC_endXYZT[idxPart]);
      momentumStart.GetXYZT(MC_startMomentum[idxPart]);
      momentumEnd.GetXYZT(MC_endMomentum[idxPart]);
      MC_truthlength[idxPart] = truthLength(particle);

      //if(isFiducial){
      h_MC_length->Fill(MC_truthlength[i]);
      // if (MC_vertex[i][0] < 1500 || MC_vertex[i][1] < 1500 || MC_vertex[i][2] < 1500 || MC_vertex[i][3] < 1500 || MC_vertex[i][4] < 1500 ) {
      //   MC_vertex[i] = MC_vertex[i]
      // } else MC_vertex[i] = -100;
      if(MC_pdg[idxPart] == -13 || MC_pdg[idxPart] == 321) {
        h_MC_length_Muons->Fill(MC_truthlength[idxPart]);
        if (MC_Endprocess[idxPart] =="Decay" || MC_Endprocess[idxPart] == "FastScintillation") {
          // for(int iDaughter=0; iDaughter < particle->NumberDaughters(); iDaughter++)
          //   {
          //    auto d_search = particleMap.find(particle->Daughter(iDaughter));
          //    auto const & daughter = *((*d_search).second);
          //    cout << "Daughters: " << daughter.PdgCode() << endl;
          //  }
          // if(MC_pdg[i] == -13 || MC_pdg[i] == 321){ 
          h_MC_PDG->Fill( MC_pdg[idxPart]);
        } else h_MC_PDGM->Fill( MC_pdg[idxPart]);
      } // if muon or kaon

      // if (MC_Endprocess[i] =="Decay" || MC_Endprocess[i] == "FastScintillation"){
      //   if (MC_pdg[i] ==  -13){
      //     // cout<<"Muon Process: " << MC_process[i] << endl;
      //     // cout<<"Muon EndProcess: " <<MC_Endprocess[i] <<endl;
      //     MC_length_Muons[i] = MC_truthlength[i];
      //     h_MC_length_Muons->Fill(MC_truthlength[i]);
      //   } else MC_length_Muons[i] = -10;
      //   h_MC_length_Muons->Fill(-10);}
      if (MC_pdg[idxPart] == 211) {
        // cout<<"Pion Process: " << MC_process[i] << endl;
        // cout<<"Pion EndProcess: " <<MC_Endprocess[i] <<endl;
        MC_length_Pions[idxPart] = MC_truthlength[idxPart];
        h_MC_length_Pions->Fill(MC_truthlength[idxPart]);
      } else MC_length_Pions[idxPart] = -10;
      //h_MC_length_Pions->Fill(-10);}
      /*
      if (MC_pdg[i] == -211){
      h_MC_length_NegPions->Fill(MC_truthlength[i]);
      }

      if (MC_pdg[i] == 321){
      //cout <<"Kaon EndProcess: " <<MC_Endprocess[i] <<endl;
      //cout <<"Kaon Process: " <<MC_process[i] << endl;
      MC_length_Kaons[i] = MC_truthlength[i];
      h_MC_length_Kaons->Fill(MC_truthlength[i]);
      }else{MC_length_Kaons[i] = -10;}
      if(MC_pdg[i] == 2112){
      //cout<<"Proton Process: " << MC_process[i] <<endl;
      //cout<<"Proton EndProcess: " <<MC_process[i] <<endl;
      h_MC_length_Protons->Fill(MC_truthlength[i]);
      }
      if(MC_pdg[i] == -11 || MC_pdg[i] == 11){
      //cout<<"Electron Process: " << MC_process[i] <<endl;
      //cout<<"Electron EndProcess: " <<MC_Endprocess[i] <<endl;
      h_MC_length_Electrons->Fill(MC_truthlength[i]);
      }
      if(MC_pdg[i] == 12 || MC_pdg[i] == 14 || MC_pdg[i] == -14 || MC_pdg[i] == 22){
      h_MC_length_Neutrinos->Fill(MC_truthlength[i]);
      }
      //h_MC_length_Kaons->Fill(-10);}
      if(MC_pdg[i] == -13 || MC_pdg[i] == 211 || MC_pdg[i] == 321 || MC_pdg[i] == 2212){
      MC_length_KMP[i] = MC_truthlength[i]; 
      }else{MC_length_KMP[i] = -10;}
      // }
      */      
      // if (particle->PdgCode() == 13 || particle->PdgCode() == 211 || particle->PdgCode() == 321 || particle->PdgCode() == 2212) {
      //   MC_pdg[i] = particle->PdgCode();
      // } else {MC_pdg[i] = 2000;}
       
      // cout<<"Process: "<<MC_process[i]<<endl;
      // cout<<"End Process: "<<MC_Endprocess[i]<<endl;
      // h_MC_end_Momentum->Fill(particle->EndPosition()(MC_endMomentum[i]));
      // h_MC_Momentum->Fill(MC_Momentum[i]);
       
      if (MC_pdg[idxPart] == 321) {
         
        h_MC_length_Kaons->Fill(MC_truthlength[idxPart]);
      
        // cout << "True PDG " << MC_pdg[i] << endl;
        // h_MC_end_Momentum_Kaons->Fill(momentumEnd.GetXYZT(MC_endMomentum[i]));
        // h_MC_Momentum_Kaons->Fill(MC_Momentum[i]);
        cout << "Kaon Track Length " << MC_truthlength[idxPart] << endl;
        // cout << "Kaon Momentum "  <<particle->Momentum().E() << endl;
        // cout << "Start Pos. X, Y, Z " << particle->Position().X() <<", " << particle->Position().Y() <<", " <<  particle->Position().Z() << endl;
        // cout << "End Pos.   X, Y, Z " << particle->EndPosition().X() <<", " << particle->EndPosition().Y() <<", " <<  particle->EndPosition().Z() << endl;
        // cout << "Kaon Process " << particle->Process() << endl;
        // cout << "Kaon End Process " << particle->EndProcess() << endl;

        for (int iDaughter=0; iDaughter < particle->NumberDaughters(); iDaughter++)
        {
          auto d_search = particleMap.find(particle->Daughter(iDaughter));
          auto const & daughter = *((*d_search).second);
          cout << "Daughters: " << daughter.PdgCode() << endl;

          if (daughter.PdgCode() == -13) {
            h_MC_DaughterMuons->Fill(MC_truthlength[idxPart]);
            h_MC_MotherKaons->Fill(MC_truthlength[idxPart]);
          }
        } // for daughter
        // Kaon_Tvertex[i] = MC_vertex[i];
      } // if kaon  

      if (MC_pdg[idxPart] == -13) {
        // h_MC_end_Momentum_Muons->Fill(momentumEnd.GetXYZT(MC_endMomentum[i]));
        // h_MC_Momentum_Muons->Fill(MC_Momentum[i]);
        // cout << "True PDG " << MC_pdg[i] << endl;
        cout << "Muon Track Length " << MC_truthlength[idxPart] << endl;
        // cout << "Muon Momentum"  <<particle->Momentum().E() <<endl;
        // cout << "Start Pos. X, Y, Z " << particle->Position().X() <<", " << particle->Position().Y() <<", " <<  particle->Position().Z() << endl;
        // cout << "End Pos.   X, Y, Z " << particle->EndPosition().X() <<", " << particle->EndPosition().Y() <<", " <<  particle->EndPosition().Z() << endl;
        // cout << "Muon Process" << particle->Process() << endl;
        // cout << "Muon End Process" << particle->EndProcess() << endl;


        // for(int iDaughter=0; iDaughter < particle->NumberDaughters(); iDaughter++)
        // {
        //  auto d_search = particleMap.find(particle->Daughter(iDaughter));
        //  auto const & daughter = *((*d_search).second);
        //  cout << "Daughters: " << daughter.PdgCode() << endl;
        // }
      } // if muon

      //}
      // Kaon_Tvertex[i] = MC_vertex[i];
      // }
      // if(MC_vertex[i] < 1500 && MC_vertex[i] > -1000) {
      //   MC_Vertex[i] = MC_vertex[i];
      // }
   
      // copy all trajectory point positions
      int tot_pts = 0;
      double* trj_xyz = (double*) MC_xyz[idxPart];
      for (auto ipt = particle->Trajectory().begin(); ipt != particle->Trajectory().end(); ipt++) {
        if (tot_pts >= MAX_MC_TRAJ_PTS) break;
        ipt->first.Vect().GetXYZ(trj_xyz+tot_pts*3);
        tot_pts++;
      }
      MC_npoints[idxPart] = tot_pts;

      /* Incomplete attempt

      //======================================================================
      // Connect true mc particle with reco tracks
      //======================================================================
      // TODO: use this to get the associated reco tracks/pfparticles
      //fTruthUtil.GetPFParticleListFromMCParticle();

      // or use particle to hits to tracks - this way have the relevant hit count already
      // need to store track ID and then do the look up in the stored tracks to store the index

      ServiceHandle<BackTrackerService> bt_svc;

      auto IDEs =  bt_svc -> TrackIdToSimIDEs_Ps(particle->TrackID());
      auto hits =  bt_svc -> TrackIdToHits_Ps(particle->TrackID());
      auto rtl = fTruthUtil.GetRecoTrackListFromMCParticle(*particle, event, fTrackModuleLabel);
      FindManyP<Hit> track_hits(trackListHandle, event, fTrackModuleLabel);

      for (auto rtrk: rtl) {

      }
      */

      //MC_Prange[i] = myPrange(MC_truthlength[i]);
      ++idxPart; //particle index

    } // for MCParticle
    MC_process.resize(idxPart);
    MC_Endprocess.resize(idxPart);

    // cout<<"dbg abc: "<<idbg++<<endl; 
    isFiducial = insideFV(MC_vertex);
    if (!isFiducial) return;


    //========================================================================
    //========================================================================
    // Reco  stuff
    //========================================================================
    //========================================================================

    // auto const &trackListHandle = event.getValidHandle<vector<Track>>(fPandoraTrackModuleLabel);
    // auto const &PandoraTrackObjs = *TrackListHandle; 
    // auto const &PfparticleHandle = event.getValidHandle<vector<recob::PFParticle>>(fPandoraLabel);
    // auto const &PandoraPfparticleObjs = *PfparticleHandle;

    // //============= ASSOCIATIONS PFP+TRACK, PFP+SHOWER, PFP+VERTEX, TRACK+CALO, TRACK+PID ==================
    // FindManyP<Track> PandoraPfp_trk(PfparticleHandle, event, fPandoraTrackModuleLabel);
    // FindManyP<recob::Shower> PandoraPfp_Swr(PfparticleHandle, event, fPandoraShowerModuleLabel);
    // FindManyP<Vertex> PandoraPfp_vtx(PfparticleHandle, event, PandoraLabel);
    // FindManyP<Calorimetry> PandoraTrk_Cal(trackListHandle, event, fPandoraCaloModuleLabel);
    // FindManyP<ParticleID> PandoraTrk_Pid(trackListHandle, event, fPandoraPIDModuleLabel);

    // //==================SAVE ONLY PARTICLE INFORMATION==========================================================
    // nPFParticle = PandoraPfparticleObjs.size();
    // cout <<"PFP Size " << nPFParticle << endl;
    // for (size_t iPfp = 0; iPfp < PandoraPfparticleObjs.size(); iPfp++)
    // {
    //   PFP_PdgCode.push_back(PandoraPfparticleObjs[iPfp].PdgCode());
    //   cout<<"PFP PDG " << PandoraPfparticleObjs[iPfp].PdgCode() << endl;
    //   if (PandoraPfparticleObjs[iPfp].IsPrimary()) 
    //   {
    //     nPrimaries++;
    //     cout<<"Is Primary " << PandoraPfparticleObjs[iPfp].IsPrimary()<< endl;
    //   }
    //   cout<<"Is Primary " << PandoraPfparticleObjs[iPfp].IsPrimary()<< endl;
    //   nDaughters.push_back(PandoraPfparticleObjs[iPfp].NumDaughters());
    //   cout<<"PFP Num. Daughters " << PandoraPfparticleObjs[iPfp].NumDaughters() <<endl;
    //   PFParticleID.push_back(PandoraPfparticleObjs[iPfp].Self());
    //   cout<<"PFP Self " <<PandoraPfparticleObjs[iPfp].Self() << endl;
    //   PFP_Parent.push_back(PandoraPfparticleObjs[iPfp].Parent());
    //   cout<<"PFP Parent " <<PandoraPfparticleObjs[iPfp].Parent() << endl;
    //   cout<"Parent PDG "  <<PandoraPfparticleObjs[PandoraPfparticleObjs[iPfp].Parent()].PdgCode() << endl; 
    // }   

    // for (size_t iPfp = 0; iPfp < PandoraPfparticleObjs.size(); iPfp++)
    // {
    //   auto const &Pfp_trk = PandoraPfp_trk.at(iPfp);

    //   PFPNumber.push_back(iPfp); // save particle flow obj ID
    //   PFPTrkSize.push_back(Pfp_trk.size());
    //   for (size_t i = 0; i < Pfp_trk.size(); i++)
    //   {
    //     auto const thisTrkPfp = Pfp_trk[i];
    //     PFPTrkLengths.push_back(thisTrkPfp->Length());
    //     PFPStartMomentum.push_back(thisTrkPfp->StartMomentum());
    //     PFPEndMomentum.push_back(thisTrkPfp->EndMomentum());
             
    //     PFPStartPoint.emplace_back(thisTrkPfp->Start().X(), thisTrkPfp->Start().Y(), thisTrkPfp->Start().Z());
    //     PFPEndPoint.emplace_back(thisTrkPfp->End().X(), thisTrkPfp->End().Y(), thisTrkPfp->End().Z());
    //     PFPStartDirection.emplace_back(thisTrkPfp->StartDirection().X(), thisTrkPfp->StartDirection().Y(), thisTrkPfp->StartDirection().Z());
                 
    //     double Length = thisTrkPfp->Length();
    //     cout << "Track Length: " << Length <<endl;
               
    //     if(PandoraPfparticleObjs[iPfp].PdgCode() == 321 || PandoraPfparticleObjs[iPfp].PdgCode() == 321) 
    //     {
    //       cout <<"PFP Kaon Length " << thisTrkPfp->Length() << endl;
    //       cout <<"Kaon Primary?? " << PandoraPfparticleObjs[iPfp].IsPrimary() << endl;
    //     }
    //     if(PandoraPfparticleObjs[iPfp].PdgCode() == 13)
    //     {
    //       cout <<"PFP Muon Length " << thisTrkPfp->Length() << endl;
    //       cout <<"PFP Muon Start X " << thisTrkPfp->Start().X() << endl;
    //       cout <<"PFP Muon Start Y " << thisTrkPfp->Start().Y() <<endl;
    //       cout <<"PFP Muon Start Z " << thisTrkPfp->Start().Z() <<endl;
    //       cout <<"PFP Muon End X " << thisTrkPfp->End().X() << endl;
    //       cout <<"PFP Muon End Y " << thisTrkPfp->End().Y() <<endl;
    //       cout <<"PFP Muon End Z " << thisTrkPfp->End().Z() <<endl;

    //       h_track_length_Muons->Fill(thisTrkPfp->Length());
    //       cout << "Muon Primary?? " << PandoraPfparticleObjs[iPfp].IsPrimary() << endl; 
    //     }

    //     //=========================ASSOCIATION OF TRACKS WITH CALO ====================================
    //     bool tmp_isContained = insideFV(thisTrkPfp->End());
    //     Track_isContained.push_back(tmp_isContained);
                      
    //     auto const &Trk_cal = PandoraTrk_Cal.at(i);
    //     if(Trk_cal.empty()) continue;
    //     nTrk_Cal.push_back(Trk_cal.size());
    //     for(size_t p = 0; p < Trk_cal.size(); p++)
    //     {
    //       auto const thisCal = Trk_cal[p];
    //       if(!thisCal->PlaneID().isValid) continue;            
    //     }
    //   }
    // }
      
    int Muon = 0;
    int Pion =0;
    // this is a track base analysis so it must be at least one track
    Handle< vector<Track> > trackListHandle;
    vector<Ptr<Track>> tracklist;
    if (event.getByLabel(fTrackModuleLabel, trackListHandle))
      fill_ptr_vector(tracklist, trackListHandle);
    cout << "here n_recoTracks is " << n_recoTracks << " and tracklist.size() is " << tracklist.size() << endl;
    n_recoTracks = tracklist.size();
    cout << "now n_recoTracks is " << n_recoTracks << " and tracklist.size() is " << tracklist.size() << endl;
    if(n_recoTracks > MAX_TRACKS || n_recoTracks == 0) {
      n_recoTracks = 0; // make sure we don't save any fake data
      return;
    }
    
    FindManyP<Vertex> trk_from_vtx(trackListHandle, event, "pmtrack");
    // FindManyP<Vertex> trk_from_vtx(trackListHandle,event,"pandora");
    // FindManyP<Vertex> trk_from_vtx(trackListHandle,event,fTrackModuleLabel);
    cout << "after getting vertex assns n_recoTracks is " << n_recoTracks << endl;

    Handle< vector<Vertex> > vtxListHandle;
    vector<Ptr<Vertex>> vtxlist;
    if(event.getByLabel("pmtrack", vtxListHandle))
    // if(event.getByLabel("pandora", vtxListHandle))
    // if(event.getByLabel(fTrackModuleLabel, vtxListHandle))
      fill_ptr_vector(vtxlist, vtxListHandle);
    cout << "after getting vertices n_recoTracks is " << n_recoTracks << endl;

    n_vertices = vtxlist.size();
    if (n_vertices != 0) {
      for (int i = 0; i < n_vertices; ++i) {
        double tmp_vtx[3] ={-999.0,-999.0,-999.0};
        vtx_ID[i] = vtxlist[i]->ID();
        vtxlist[i]->XYZ(tmp_vtx);
        for (int j = 0; j < 3; ++j) vertex[i][j] = tmp_vtx[j];
      }
    }

    cout << "before getting hit assns n_recoTracks is " << n_recoTracks << endl;
    FindManyP<Hit> track_hits(trackListHandle, event, fTrackModuleLabel);
    art::FindMany<Calorimetry> reco_cal(trackListHandle, event, fCaloModuleLabel);
    trkf::TrackMomentumCalculator trackP;
    cout << "after getting hit assns n_recoTracks is " << n_recoTracks << endl;

    art::FindMany<ParticleID> reco_PID(trackListHandle, event, fPidModuleLabel);
    cout << "after getting particle assns n_recoTracks is " << n_recoTracks << endl;

    Handle<vector<Hit>> HitHandle;
    vector<Ptr<Hit>> all_hits;
    if(event.getByLabel(fHitModuleLabel, HitHandle))
      fill_ptr_vector(all_hits, HitHandle);
    cout << "after getting hits n_recoTracks is " << n_recoTracks << endl;

    for (int i = 0; i < n_recoTracks; ++i) {
      // if (n_recoTracks != 2) { cout << "More Than 2 Tracks!! " << endl; continue;}
      Ptr<Track> track = tracklist[i];

      // vtx associations
      cout << "size of tracklist is " << tracklist.size() << " n_recoTracks is " << n_recoTracks
        << ", size of trk_from_vtx is " << trk_from_vtx.size() << ", i is " << i << endl;
      vector<Ptr<Vertex>> vtxs = trk_from_vtx.at(i);
      for (size_t j = 0; j < vtxs.size(); ++j) {
        Ptr<Vertex> vtx = vtxs[j];
        vtxID_trk[i][j] = vtx->ID();
      }
      track_length[i] = track->Length();

      const TVector3 tmp_track_vtx(track->Vertex().x(),
        track->Vertex().y(),
        track->Vertex().z());
      const TVector3 tmp_track_end(track->End().x(),
        track->End().y(),
        track->End().z());
      track_vtx[i][0] =tmp_track_vtx[0];
      track_vtx[i][1] =tmp_track_vtx[1];
      track_vtx[i][2] =tmp_track_vtx[2];
      track_vtx[i][3] = -999.0;

      track_end[i][0] =tmp_track_end[0];
      track_end[i][1] =tmp_track_end[1];
      track_end[i][2] =tmp_track_end[2];
      track_end[i][3] = -999.0;

      const TVector3 tmp_vtx_dir(track->VertexDirection().x(),
        track->VertexDirection().y(),
        track->VertexDirection().Z());
      track_vtxDir[i][0] = tmp_vtx_dir[0];
      track_vtxDir[i][1] = tmp_vtx_dir[1];
      track_vtxDir[i][2] = tmp_vtx_dir[2];

      track_ID[i] = track->ID();
  
      // changed 13 to 20 below
      track_Prange[i] = trackP.GetTrackMomentum(track_length[i],100);
      double trk_end[4] ={tmp_track_end[0],tmp_track_end[1],tmp_track_end[2],-999};
      bool track_isInside = insideFV( trk_end );
      // check if the track ends within the FV
      if (track_isInside) track_isContained[i] = 1;
      else track_isContained[i] =0;
      //calculate PID
      cout << "size of Reco PID is " << reco_PID.size() << ", n_RecoTracks is " << n_recoTracks
        << ", size of Reco Cal is " << reco_cal.size() << endl;
      vector<const Calorimetry*> trk_cal = reco_cal.at(i);
      vector<const ParticleID*> trk_pid = reco_PID.at(i);
      int plane0 = trk_pid[0]->Ndf();
      int plane1 = trk_pid[1]->Ndf();
      int plane2 = trk_pid[2]->Ndf();
      int best_plane = -1;
      int most_ndf = max({plane0, plane1, plane2});
      for( size_t p =0; p<3; p++) if( most_ndf == trk_pid[p]->Ndf() ) best_plane = p;
      // cout << "BEST PLANE"  << best_plane <<endl;
      track_PID_pdg[i][0] = trk_pid[0]->Pdg();
      track_PID_pdg[i][1] = trk_pid[1]->Pdg();
      track_PID_pdg[i][2] = trk_pid[2]->Pdg();

      if(track_PID_pdg[i][0] == 321 || track_PID_pdg[i][0] == -13 || track_PID_pdg[i][0] == 13 || track_PID_pdg[i][0] == 11 || track_PID_pdg[i][0] == 211){
        track_PID_pdg[i][0] = track_PID_pdg[i][0];
      }else{track_PID_pdg[i][0] = -23; }
      if(track_PID_pdg[i][1] == 321 || track_PID_pdg[i][1] == -13 || track_PID_pdg[i][1] == 13 || track_PID_pdg[i][0] == 11 || track_PID_pdg[i][1]){
        track_PID_pdg[i][1] = track_PID_pdg[i][1];
      }else{track_PID_pdg[i][1] = -23; }
      if(track_PID_pdg[i][2] == 321 || track_PID_pdg[i][2] == -13 || track_PID_pdg[i][2] == 13 || track_PID_pdg[i][0] == 11 || track_PID_pdg[i][2]){
        track_PID_pdg[i][2] = track_PID_pdg[i][2];
      }else{track_PID_pdg[i][2] = -23; }

      // if (track_PID_pdg[i][0]<30){
      //   track_PID_pdg[i][0] = track_PID_pdg[i][0];
      // } else { track_PID_pdg[i][0] = 28; }
      // if (track_PID_pdg[i][1]<30){       
      //   track_PID_pdg[i][1] = track_PID_pdg[i][1];
      // } else { track_PID_pdg[i][1] = 28; }      
      // if (track_PID_pdg[i][2]<30){
      //   track_PID_pdg[i][2] = track_PID_pdg[i][2];
      // } else { track_PID_pdg[i][2] = 28; }

      // if (n_recoTracks !=2) continue; 
      h_track_length->Fill(track->Length());
      cout << "Reco Track Lengths " << track->Length() << endl;
   
      if (track_PID_pdg[i][0] == 321 || track_PID_pdg[i][1] == 321  || track_PID_pdg[i][2] == 321  ){ 
      // if (MC_Endprocess[i] =="Decay" || MC_Endprocess[i] == "FastScintillation"){
        // cout <<"Reco PID " << track_PID_pdg[i][0]  <<endl;
        // cout <<"Reco PID " << track_PID_pdg[i][1]  <<endl;
        // cout <<"Reco PID " << track_PID_pdg[i][2]  <<endl;
        // cout<<"Process Reco"<<MC_process <<endl;
        // cout <<"Reco PID 321" <<endl;
        cout <<"Reco Kaon Track Length " << track_length[i] << endl;
        // cout <<"Reco Vertex " << track_vtx[i][1] << endl;
        // cout<<"Kaon Reco Track Vertex " << track->Vertex().x() << ", " << track->Vertex().y() << ", " << track->Vertex().z() << endl;
        // cout <<"Kaon Reco Track End "   << track->End().x() << ", " << track->End().y() << ", " << track->End().z() << endl; 

        if (-352.662 > track->End().x() || 352.662 < track->End().x()|| -593.924 > track->End().y() || 593.924 < track->End().y() || 9.12375 > track->End().z() || 1383.924 < track->End().z())
          h_track_length_Walls_Kaons->Fill(track->Length());
        else
          h_track_length_NonWalls_Kaons->Fill(track->Length());

        // cout <<"Reco Decay Vertex " << decayVtx[i][3] << endl;
        // cout <<"Reco DE/dx " << track_dE_dx[i] << endl;
        // cout <<"Vertex Difference" <<MC_vertex[i] - track_vtx[i][1] << endl;
        // Kaon_Rvertex[i] = track_vtx[i][1];
           
        if (track->Length() > 0.0) {
          track_length_Kaons[i] = track->Length();
          h_track_length_Kaons->Fill(track->Length()); 
        } else {
          track_length_Kaons[i] = -10;
          //h_track_length_Kaons->Fill(-10);
        }
      } // if reco kaon

      // if (track_PID_pdg[i][0] == 13 || track_PID_pdg[i][1] == 13 || track_PID_pdg[i][2] == 13){
      // if(MC_Endprocess[i] =="Decay" || MC_Endprocess[i] == "FastScintillation"){
      if (track_PID_pdg[i][0] == 13 || track_PID_pdg[i][1] == 13  || track_PID_pdg[i][2] == 13 ) {
        Muon++;
        cout << "Muon Count " << Muon << endl;
        if (Muon > 1) {h_track_length_Kaons->Fill(track->Length()); continue;}
        track_length_Muons[i] = track->Length(); 
        if (track_length_Muons[i] == track_length_Kaons[i]) {cout << "Double Count!!" << endl; continue;}
        h_track_length_Muons->Fill(track->Length());
        cout << "Reco Muon Track Length " << track_length[i] << endl;
        // cout<<" Muon Reco Track Vertex " << track->Vertex().x() << ", " << track->Vertex().y() << ", " << track->Vertex().z() << endl;
        // cout <<"Muon Reco Track End "   << track->End().x() << ", " << track->End().y() << ", " << track->End().z() << endl;

        if (-352.662 > track->End().x() || 352.662 < track->End().x()|| -593.924 > track->End().y() || 593.924 < track->End().y() || 9.12375 > track->End().z() || 1383.924 < track->End().z())
          h_track_length_Walls_Muons->Fill(track->Length());
        else
          h_track_length_NonWalls_Muons->Fill(track->Length());

        if (track->Length() < 40) { 
          h_shortMuonsX->Fill(track->End().x());
          h_shortMuonsY->Fill(track->End().y());
          h_shortMuonsZ->Fill(track->End().z());
        } else {
          track_length_Muons[i] = -10;
          // h_track_length_Muons->Fill(-10);
        }
      } // if reco muon

      if (track_PID_pdg[i][0] == 211 || track_PID_pdg[i][1] == 211  || track_PID_pdg[i][2] == 211) {
        track_length_Pions[i] = track->Length();
        Pion++;
        cout << "Pion Count " << Pion << endl;
        // cout<<"Reco Pion Track Length " << track_length[i] << endl;
        if (track_length_Pions[i] == track_length_Kaons[i]) {cout << "Double Count!!" << endl; continue;}
        if (track_length_Pions[i] == track_length_Muons[i]) {cout << "Double Count!!" << endl; continue;}
        if (track_length_Pions[i] > 45) {h_track_length_Muons->Fill(track->Length()); continue;}
        cout << "Reco Pion Track Length " << track_length[i] << endl;
        h_track_length_Pions->Fill(track->Length());
      } else {
        track_length_Pions[i] = -10;
        //h_track_length_Pions->Fill(-10);
      }

      if (track_PID_pdg[i][0] == 2212 || track_PID_pdg[i][1] == 2212 || track_PID_pdg[i][2] || track_PID_pdg[i][0] == 13 || track_PID_pdg[i][1] == 13 || track_PID_pdg[i][2] == 13 || track_PID_pdg[i][0] == 211 || track_PID_pdg[i][1] == 211 || track_PID_pdg[i][2] == 211 || track_PID_pdg[i][0] == 211 || track_PID_pdg[i][1] == 211 || track_PID_pdg[i][2] == 211)
        track_length_KMP[i] = track->Length();
      else
        track_length_KMP[i] = -10;

      // cout << "Reco PID 1" << track_PID_pdg[i][1] << endl;
      // cout << "Reco PID 2" << track_PID_pdg[i][2] << endl;
      // cout << "Reco Track Length" << track_length[i] << endl;

      track_KE[i][0] = trk_cal[0]->KineticEnergy();
      track_KE[i][1] = trk_cal[1]->KineticEnergy();
      track_KE[i][2] = trk_cal[2]->KineticEnergy();

      vector<double> PIDAval;
      vector<double> chi2;
      PIDAcal(trk_cal, PIDAval);
      int idx = 0;
      for (auto const& val : PIDAval) {
        track_PIDA[i][idx] = val;
        idx++;
      }
      track_bestplane[i] = best_plane;

      // save dE/dx & dQ/dx
      for (unsigned iplane = 0; iplane < 3; iplane++) {
        int npts = trk_cal[iplane]->dEdx().size();
        if (npts > MAX_CALO_PTS)
        npts = MAX_CALO_PTS;
        if (best_plane >= 0 && iplane == (unsigned)best_plane) {
          n_cal_points[i] = npts;
          // if (npts < 0) n_cal_points[i] = 0;
        }

        n_cal_points_byplane[i][iplane] = npts;

        const vector<float> &dqdx = trk_cal[iplane]->dQdx();
        const vector<float> &dedx = trk_cal[iplane]->dEdx();
        const vector<float> &resr = trk_cal[iplane]->ResidualRange();
        const vector<float> &pitch = trk_cal[iplane]->TrkPitchVec();
        auto &xyz = trk_cal[iplane]->XYZ();

        if (best_plane >= 0 && iplane == (unsigned)best_plane) {
          copy(dqdx.begin(), dqdx.end(), track_dQ_dx[i]);
          copy(dedx.begin(), dedx.end(), track_dE_dx[i]);
          copy(resr.begin(), resr.end(), track_range[i]);
          copy(pitch.begin(), pitch.end(), track_pitch[i]);
        }
        copy(dqdx.begin(), dqdx.end(), track_dQ_dx_byplane[i][iplane]);
        copy(dedx.begin(), dedx.end(), track_dE_dx_byplane[i][iplane]);
        copy(resr.begin(), resr.end(), track_range_byplane[i][iplane]);
        copy(pitch.begin(), pitch.end(), track_pitch_byplane[i][iplane]);

        // save calo point's XYZ coords
        double* coords = (double*)track_calo_xyz_byplane[i][iplane];
        for(int j = 0; j < npts; j++) {
          (coords+j*3)[0] = xyz[j].X();
          (coords+j*3)[1] = xyz[j].Y();
          (coords+j*3)[2] = xyz[j].Z();
        }
      } // for plane

      //truth matcher
      double tmpEfrac = 0;
      auto const clockData = ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
      const MCParticle *particle;
      double tmpComplet = 0;
      vector<Ptr<Hit>> all_trackHits = track_hits.at(i);
      truthMatcher(clockData, all_hits, all_trackHits, particle, tmpEfrac, tmpComplet);

      if (!particle) {
        cout << "matching was unable to find a particle" << endl;
        continue;
      }
      track_mcID[i] = particle->TrackId();
      track_mcPDG[i] = particle->PdgCode();
      track_Efrac[i] = tmpEfrac;
      track_complet[i] = tmpComplet;
    
    } // Reco Track Loop

    //cout<<"dbg abcd: "<<idbg++<<endl;

    auto const clockData = ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    truthMatcher2(clockData,tracklist, track_hits, mc_id2idx);

    // CNN dacayID vertex
    Handle<vector<Vertex>> dcy_vtxHandle;
    vector<Ptr<Vertex>> dcy_vtxlist;
    if(event.getByLabel(fPointIDModuleLabel,dcy_vtxHandle))
      fill_ptr_vector(dcy_vtxlist, dcy_vtxHandle);
    n_decayVtx= dcy_vtxlist.size();

    //FindManyP<Track> decay_tracklist(dcy_vtxHandle, event, fPointIDModuleLabel);
    if (n_decayVtx !=0) {
      for (int i=0; i< n_decayVtx; ++i) {
        double tmp_vtx[3] = {-999.0,-999.0,-999.0};
        dcy_vtxlist[i]->XYZ(tmp_vtx);
        for (int j=0; j<3; ++j) decayVtx[i][j] = tmp_vtx[j];
        // vector<Ptr<Track>>  decay_track = decay_tracklist.at(i);
        // cout<<"how many tracks? "<<decay_track.size()<<endl;
      }
    }

    // CNN Em-trk hits
    // Imported code from PointIdEffTest_module.cc
    // to included hit and CNN output info into the analysis module
    // to quantify how much EM activity we have in order to reduce background

    /**
     * Needed to comment the followin out in order to be able to run
     * with latest reconstruction. Exception was thrown complaining
     * about different sizes of feature vector and data products.
     *
     * When rerunning linecluster algorithm and not blurredcluster,
     * emtrkmichelid, and emshower. This causes mismatch between
     * associations for TMVA related products which are tight to
     * linecluster products.
     **/
    /* 
    anab::MVAReader<Hit, MVA_LENGTH> hitResults(event, fNNetModuleLabel);                     // hit-by-hit outpus just to be dumped to file for debugging
    auto cluResults = anab::MVAReader<recob::Cluster, MVA_LENGTH>::create(event, fNNetModuleLabel);  // outputs for clusters recovered in not-throwing way
    int trk_idx = hitResults.getIndex("track");
    int em_idx = hitResults.getIndex("em");
    //int michel_idx = hitResults.getIndex("michel");
    Em_ch =0.0;
    Em_e =0.0;
    trk_e =0.0;
    if (cluResults) {
      const FindManyP<Hit> hitsFromClusters(cluResults->dataHandle(), event, cluResults->dataTag());
      for (size_t c = 0; c < cluResults->size(); ++c) {
        const recob::Cluster & clu = cluResults->item(c);
        if(clu.Plane().Plane != fView) continue;
        const vector< Ptr<Hit> > & hits = hitsFromClusters.at(c);
        vector< anab::FeatureVector<MVA_LENGTH> >  hit_outs = hitResults.outputs();
        //EM hits
        for (auto const & h : hits) {
          auto const & vout = hit_outs[h.key()];
          double p_trk_or_sh = vout[trk_idx] + vout[em_idx];
          if (p_trk_or_sh > 0) {
            double PidValue = vout[trk_idx] / p_trk_or_sh;
            if (PidValue < fPidValue) {
              Em_ch += h->SummedADC()* fCalorimetryAlg.LifetimeCorrection( h->PeakTime() );
              Em_e  += (fCalorimetryAlg.ElectronsFromADCArea( h->Integral(), h->WireID().Plane) * fCalorimetryAlg.LifetimeCorrection( h->PeakTime() ) ) / util::kGeVToElectrons;
              //Michel hits
              //if( vout[michel_idx] > 0.1 ) //temporay
              //Emichel_e  += fCalorimetryAlg.ElectronsFromADCArea( h->Integral(), h->WireID().Plane) * fCalorimetryAlg.LifetimeCorrection( h->PeakTime() );
            }
            else  trk_e += (fCalorimetryAlg.ElectronsFromADCArea( h->Integral(), h->WireID().Plane) * fCalorimetryAlg.LifetimeCorrection( h->PeakTime() ) ) / util::kGeVToElectrons;
          }
        }
      }
    }
     
    */
    /**
     * End of temporary commented out code.
     **/

    // Showers... for background rejection
    Handle<vector<recob::Shower>> showerHandle;
    vector<Ptr<recob::Shower>> showerlist;
    if(event.getByLabel(fShowerModuleLabel,showerHandle))
      fill_ptr_vector(showerlist, showerHandle);
    n_recoShowers= showerlist.size();
    if( n_recoShowers != 0 ) {
      for (int i=0; i<n_recoShowers && i< MAX_SHOWERS; ++i) {
        Ptr<recob::Shower> shower = showerlist[i];
        sh_direction_X[i] = shower->Direction().X();
        sh_direction_Y[i] = shower->Direction().Y();
        sh_direction_Z[i] = shower->Direction().Z();
        sh_start_X[i] = shower->ShowerStart().X();
        sh_start_Y[i] = shower->ShowerStart().Y();
        sh_start_Z[i] = shower->ShowerStart().Z();
        sh_bestplane[i] = shower->best_plane();
        sh_length[i] = shower->Length();
        for( size_t j =0; j<shower->Energy().size(); j ++) sh_energy[i][j] = shower->Energy()[j];
        for( size_t j =0; j<shower->MIPEnergy().size(); j++) sh_MIPenergy[i][j] = shower->MIPEnergy()[j];
        for( size_t j =0; j<shower->dEdx().size(); j++) sh_dEdx[i][j] = shower->dEdx()[j];
      }
    }

    if (fSaveHits) {
      //Hits
      n_recoHits= all_hits.size();
      if (n_recoHits != 0)
      for (int i = 0; i < n_recoHits && i < MAX_HITS ; ++i) {//loop over hits
        hit_channel[i] = all_hits[i]->Channel();
        hit_tpc[i]   = all_hits[i]->WireID().TPC;
        hit_plane[i]   = all_hits[i]->WireID().Plane;
        hit_wire[i]    = all_hits[i]->WireID().Wire;
        hit_peakT[i]   = all_hits[i]->PeakTime();
        hit_charge[i]  = all_hits[i]->Integral();
        hit_ph[i]  = all_hits[i]->PeakAmplitude();
        hit_startT[i] = all_hits[i]->PeakTimeMinusRMS();
        hit_endT[i] = all_hits[i]->PeakTimePlusRMS();
        hit_rms[i] = all_hits[i]->RMS();
        hit_electrons[i]  = fCalorimetryAlg.ElectronsFromADCArea( all_hits[i]->Integral(), all_hits[i]->WireID().Plane);
        // fCalorimetryAlg.LifetimeCorrection(clockData, detProp, theHit./gser.TimeToCm() );
      } // for hit
    } // if saving hits
  } // function NDKAna::Process

  //========================================================================
  void NDKAna::PIDAcal(vector<const Calorimetry*> cal, vector<double> &PIDA) {
    vector<double> pida_vec_v0;
    vector<double> pida_vec_v1;
    vector<double> pida_vec_v2;

    int used_points[3] = {0,0,0};
    double tmp_pida[3];
    for (unsigned j =0; j<3; ++j) {
      for (unsigned i =0; i<cal[j]->dEdx().size(); ++i) { // loop through hits on each plane
        if (cal[j]->ResidualRange()[i] < fPIDA_endPoint) { // Only want PIDA for last x cm
          tmp_pida[j] = cal[j]->dEdx()[i]* pow(cal[j]->ResidualRange()[i], fExponentConstant );
          if (fMinPIDAValue > tmp_pida[j] || tmp_pida[j] > fMaxPIDAValue) continue;
          if (j == 0) pida_vec_v0.push_back(tmp_pida[j]);
          if (j == 1) pida_vec_v1.push_back(tmp_pida[j]);
          if (j == 2) pida_vec_v2.push_back(tmp_pida[j]);
          used_points[j]++;
        } // If ResRange < x cm
      }// Loop over hits on each plane
    }

    // for each view calculate PIDA median value
    sort(pida_vec_v0.begin(), pida_vec_v0.end());
    int size_v0 = pida_vec_v0.size();
    double median_v0 = -999;
    if (size_v0 > 0) median_v0 = size_v0 % 2 ? pida_vec_v0[size_v0 / 2] : (pida_vec_v0[size_v0 / 2 - 1] + pida_vec_v0[size_v0 / 2]) / 2;

    sort(pida_vec_v1.begin(), pida_vec_v1.end() );
    int size_v1 = pida_vec_v1.size();
    double median_v1 = -999;
    if (size_v1 > 0) median_v1 = size_v1 % 2 ? pida_vec_v1[size_v1 / 2] : (pida_vec_v1[size_v1 / 2 - 1] + pida_vec_v1[size_v1 / 2]) / 2;

    sort(pida_vec_v2.begin(), pida_vec_v2.end() );
    int size_v2 = pida_vec_v2.size();
    double median_v2 = -999;
    if (size_v2 > 0) median_v2 = size_v2 % 2 ? pida_vec_v2[size_v2 / 2] : (pida_vec_v2[size_v2 / 2 - 1] + pida_vec_v2[size_v2 / 2]) / 2;

    PIDA.push_back(median_v0);
    PIDA.push_back(median_v1);
    PIDA.push_back(median_v2);

  } // function NDKAna::PIDAcal

  //========================================================================
  void NDKAna::truthMatcher(DetectorClocksData const& clockData, 
                            vector<Ptr<Hit>> all_hits,
                            vector<Ptr<Hit>> track_hits, 
                            const MCParticle*& MCparticle, 
                            double &Efrac, 
                            double &Ecomplet)
  {
    // get the usefull tools
    ServiceHandle<BackTrackerService const> bt;
    ServiceHandle<ParticleInventoryService> part_inv;
    cout << "KANGAROO FRIENDS" <<endl;
    // fill the map of deposited energy to true track id
    map<int,double> trkID_E;
    for (size_t j = 0; j < track_hits.size(); ++j) {
      Ptr<Hit> hit = track_hits[j];
      // vector<TrackIDE> TrackIDs = bt->HitToTrackIDEs(hit);
      vector<TrackIDE> TrackIDs = bt->HitToTrackIDEs(clockData,hit);
      // cout << "TrackIDs Size " <<TrackIDs.size() << endl;
      for(size_t k = 0; k < TrackIDs.size(); k++){
        // note here that stored energy will be larger than actual
        // deposited energy, since 1 IDE will produce hits in
        // multiple wire planes!
        trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
        // cout << "Matching PID " << part_inv->TrackIdToParticle(k.trackID)<< endl;
        // MCParticle p = part_inv->TrackIdToParticle(k);
        // cout << "PdgCode Butts " << p.PdgCode() <<  endl;
      }
    }

    // total energy budget accounting
    double E_em =0.0;
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla
    //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition
    //!since we are looking for muons/pions/protons this should be enough
    if( !trkID_E.size() ) {
      MCparticle = 0;
      return; //Ghost track???
    }
    // loop over stored true track depositions
    for (map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii) {
      total_E += ii->second;
      // pick the true track with most energy deposited
      if((ii->second)>max_E){
        partial_E = ii->second;
        max_E = ii->second;
        TrackID = ii->first;
        if( TrackID < 0 ) E_em += ii->second;
      }
    }
     
    MCparticle = part_inv->TrackIdToParticle_P(TrackID);
    cout << "Track ID " << TrackID << endl;
    cout << "MCparticle Matched" <<MCparticle << endl;
      
    //In the current simulation, we do not save EM Shower daughters in GEANT. But we do save the energy deposition in TrackIDEs. If the energy deposition is from a particle that is the daughter of
    //an EM particle, the negative of the parent track ID is saved in TrackIDE for the daughter particle
    //we don't want to track gammas or any other EM activity
    if (TrackID < 0) return;
    cout <<"TrackID < 0 Cut " << endl;
    // Efrac = (partial_E+E_em)/total_E;
    Efrac = partial_E / total_E;

    // completeness
    double totenergy =0;
    for (size_t k = 0; k < all_hits.size(); ++k) {
      Ptr<Hit> hit = all_hits[k];
      vector<TrackIDE> TrackIDs = bt->HitToTrackIDEs(clockData,hit);
      for (size_t l = 0; l < TrackIDs.size(); ++l) {
        if (TrackIDs[l].trackID==TrackID) totenergy += TrackIDs[l].energy;
      }
    } // for hit

    Ecomplet = partial_E / totenergy;

  } // function NDKAna::truthMatcher

  //========================================================================
  void NDKAna::truthMatcher2(DetectorClocksData const& clockData,
                             vector<Ptr<Track> > &tracklist,
                             FindManyP<Hit> &track_hits,
                             map<int, int> &mc_id2idx)
  {
    if (!tracklist.size()) return;

    //int idbg = 0;
    //cout<<"dbg abcde: "<<idbg++<<endl;
    // get the usefull tools
    ServiceHandle<BackTrackerService> bt;
    ServiceHandle<ParticleInventoryService> part_inv;

    cout << "TruthMatcher2 Pre Mapping" << endl;
    map<int, map<int, int> > truetrk2hits; // map true track to number of hits by reco tracks
    //map<int, vector<int> > truetrk2recotrack; // map true track to ID of reco tracks

    cout <<"truthMatcher2 Inside" <<endl;
    // filling the map relating true tracks, reco tracjks and number of hits
    int size = track_hits.size();
    for (int i = 0; i < size; i++) {
      for (auto hit : track_hits.at(i)) {
        vector<TrackIDE> trkIDEs = bt->HitToTrackIDEs(clockData,hit);
        set<int> trkIDs;
        for (auto ide: trkIDEs)
          trkIDs.insert(ide.trackID); // which true tracks contributed to this hit
        for (auto trkID : trkIDs) {
          if (!truetrk2hits.count(trkID) ||
            !truetrk2hits[trkID].count(tracklist[i]->ID())) {
            truetrk2hits[trkID][tracklist[i]->ID()] = 0;
          }
          truetrk2hits[mc_id2idx[trkID]][tracklist[i]->ID()]++;
        } // for track id
      } // for hit
    } // for track
    //cout<<"dbg abcdef: "<<idbg++<<endl;

    // fill the true-reco associations
    /*
    for (auto trueit = truetrk2hits.begin(); trueit != truetrk2hits.end(); trueit++) {
      int i = 0;
      for (auto trk2nhit: trueit->second) {
        MC_reco_track_ids[trueit->first][i] = trk2nhit.first;
        MC_nhits_in_track[trueit->first][i] = trk2nhit.second;
        i++;
      }
      MC_nrecotracks[trueit->first] = i;
    }
    cout<<"dbg abcdefg: "<<idbg++<<endl;
    */
    //Issue Here

  } // function NDKAna::truthMatcher2

  //========================================================================
  double NDKAna::truthLength(const MCParticle *MCparticle) {
    //calculate the truth length considering only the part that is inside the TPC
    //Base on a peace of code from dune/TrackingAna/TrackingEfficiency_module.cc

    if(!MCparticle) return -999.0;
    double TPCLength = 0.0;
    int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
    vector<double> TPCLengthHits(numberTrajectoryPoints, 0);
    int FirstHit=0, LastHit=0;
    bool BeenInVolume = false;

    for (int MCHit = 0; MCHit < numberTrajectoryPoints; ++MCHit) {
      const TLorentzVector& tmpPosition= MCparticle->Position(MCHit);
      double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
      if (MCHit!=0) TPCLengthHits[MCHit] = sqrt( pow( (MCparticle->Vx(MCHit-1)-MCparticle->Vx(MCHit)),2)+ pow( (MCparticle->Vy(MCHit-1)-MCparticle->Vy(MCHit)),2)+ pow( (MCparticle->Vz(MCHit-1)-MCparticle->Vz(MCHit)),2));

      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
      if (tpcid.isValid) {
        // -- Check if hit is within drift window...
        geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat);
        geo::TPCGeo      const& tpc  = cryo.TPC(tpcid.TPC);
        double XPlanePosition      = tpc.PlaneLocation(0)[0];
        double DriftTimeCorrection = fabs( tmpPosition[0] - XPlanePosition ) / XDriftVelocity;
        double TimeAtPlane         = MCparticle->T() + DriftTimeCorrection;
        // if( TimeAtPlane < detprop->TriggerOffset() || TimeAtPlane > detprop->TriggerOffset() + WindowSize ) continue;
        if( TimeAtPlane < detprop.TimeOffsetZ() || TimeAtPlane > detprop.TimeOffsetZ() + WindowSize ) continue;
        LastHit = MCHit;
        if(!BeenInVolume) {
          BeenInVolume = true;
          FirstHit = MCHit;
        }
      }
    }
    for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) TPCLength += TPCLengthHits[Hit];
    return TPCLength;
  } // function NDKAna::truthMatcher2

  //========================================================================
  bool NDKAna::insideFV(double vertex[4]) {

    double x = vertex[0];
    double y = vertex[1];
    double z = vertex[2];

    if (x>fFidVolXmin && x<fFidVolXmax&&
      y>fFidVolYmin && y<fFidVolYmax&&
      z>fFidVolZmin && z<fFidVolZmax)
      return true;
    else
      return false;
  } // function NDKAna::insideFV

  //========================================================================
  void NDKAna::reset() {

    //cout<<"NDKAna::reset()"<<endl;
    Event = 0;
    Run = 0;
    SubRun = 0;

    //MC truth
    MC_Ev = 0;
    MC_nuPDG = 0;
    MC_Q2 = 0;
    MC_hit_nucleon = 0;
    MC_cc = 0;
    MCgenie_npart = 0;
    memset( MCgenie_id, 0, sizeof(MCgenie_id) );
    memset( MCgenie_pdg, 0, sizeof(MCgenie_pdg) );
    memset( MCgenie_mother, 0, sizeof(MCgenie_mother) );
    memset( MCgenie_statusCode, 0, sizeof(MCgenie_statusCode) );
    memset( MCgenie_fate, 0, sizeof(MCgenie_fate) );
    memset( MCgenie_startMomentum, 0, sizeof(MCgenie_startMomentum) );
    memset( MCgenie_endMomentum, 0, sizeof(MCgenie_endMomentum) );
    //cout<<"NDKAna::reset(): cleared MCgenie data holders"<<endl;

    memset( MC_vertex, 0, sizeof(MC_vertex) );
    //memset( MC_Vertex, 0, sizeof(MC_vertex) );
    //cout<<"NDKAna::reset(): cleared MC_vertex"<<endl;
    MC_npart = 0;
    memset( MC_id, 0, sizeof(MC_id) );
    //cout<<"NDKAna::reset(): cleared MC_id"<<endl;
    memset( MC_pdg, 0, sizeof(MC_pdg) );
    //cout<<"NDKAna::reset(): cleared MC_pdg"<<endl;
    memset( MC_mother, 0, sizeof(MC_mother) );
    memset( MC_startXYZT, 0, sizeof(MC_startXYZT) );
    memset( MC_endXYZT, 0, sizeof(MC_endXYZT) );
    memset( MC_startMomentum, 0, sizeof(MC_startMomentum) );
    memset( MC_endMomentum, 0, sizeof(MC_endMomentum) );
    memset( MC_truthlength, 0, sizeof(MC_truthlength) );
    memset( MC_length_Kaons, 0, sizeof(MC_truthlength) );
    memset( MC_length_Muons, 0, sizeof(MC_truthlength) );
    memset( MC_length_Pions, 0, sizeof(MC_truthlength) );
    memset( MC_length_KMP,  0, size(MC_truthlength) );



    memset( MC_Prange, 0, sizeof(MC_Prange) );
    memset( MC_statusCode, 0, sizeof(MC_statusCode) );
    cout <<MC_process.size()<<endl;
    //cout<<"NDKAna::reset(): to clear MC_process"<<endl;
    MC_process.clear();
    //cout<<"NDKAna::reset(): to clear MC_Endprocess"<<endl;
    MC_Endprocess.clear();
    //cout<<"NDKAna::reset(): cleared MC truth data holders"<<endl;

    //cout<<"NDKAna::reset(): about to clear new MC data holders."<<endl;
    memset( MC_nrecotracks, 0, sizeof(MC_nrecotracks) );
    memset( MC_reco_track_ids, 0, sizeof(MC_reco_track_ids) );
    memset( MC_nhits_in_track, 0, sizeof(MC_nhits_in_track) );
    memset( MC_energy_in_track, 0, sizeof(MC_energy_in_track) );
    //cout<<"NDKAna::reset(): done clearing new MC data holders."<<endl;

    n_vertices = 0;
    memset( vertex, 0, sizeof(vertex) );
    memset( vtx_ID, 0, sizeof(vtx_ID) );
    n_decayVtx = 0;
    memset( decayVtx, 0, sizeof(decayVtx) );
    n_recoTracks = 0;
    memset( track_isContained, 0, sizeof(track_isContained) );
    memset( track_ID, 0, sizeof(track_ID) );
    memset( vtxID_trk, 0, sizeof(vtxID_trk) );
    memset( track_vtxDir, 0, sizeof(track_vtxDir) );
    memset( track_vtx, 0, sizeof(track_vtx) );
    memset( track_end, 0, sizeof(track_end) );
    memset( track_length, 0, sizeof(track_length) );
    memset( track_length_Kaons, 0, sizeof(track_length_Kaons) );  
    memset( track_length_Muons, 0, sizeof(track_length_Muons) );
    memset( track_length_Pions, 0, sizeof(track_length_Pions) );
    memset( track_length_KMP,   0, sizeof(track_length_KMP)   );   

    memset( track_dir_vtx, 0, sizeof(track_dir_vtx) );
    memset( track_PIDA, 0, sizeof(track_PIDA) );
    memset( track_PID_pdg, 0, sizeof(track_PID_pdg) );
    memset( track_KE, 0, sizeof(track_KE) );
    memset( track_bestplane, 0, sizeof(track_bestplane) );
    memset( track_Prange, 0, sizeof(track_Prange) );
    memset( track_Efrac, 0, sizeof(track_Efrac) );
    memset( track_complet, 0, sizeof(track_complet) );
    memset( track_mcID, 0, sizeof(track_mcID) );
    memset( track_mcPDG, 0, sizeof(track_mcPDG) );
    memset( n_cal_points, 0, sizeof(n_cal_points) );
    memset( track_dQ_dx, 0, sizeof(track_dQ_dx) );
    memset( track_dE_dx, 0, sizeof(track_dE_dx) );
    memset( track_range, 0, sizeof(track_range) );
    memset( track_pitch, 0, sizeof(track_pitch) );
    memset( n_cal_points_byplane, 0, sizeof(n_cal_points_byplane) );
    memset( track_dQ_dx_byplane, 0, sizeof(track_dQ_dx_byplane) );
    memset( track_dE_dx_byplane, 0, sizeof(track_dE_dx_byplane) );
    memset( track_range_byplane, 0, sizeof(track_range_byplane) );
    memset( track_pitch_byplane, 0, sizeof(track_pitch_byplane) );

    n_recoHits = 0;
    memset( hit_channel, 0, sizeof(hit_channel) );
    memset( hit_tpc, 0, sizeof(hit_tpc) );
    memset( hit_plane, 0, sizeof(hit_plane) );
    memset( hit_wire, 0, sizeof(hit_wire) );
    memset( hit_peakT, 0, sizeof(hit_peakT) );
    memset( hit_charge, 0, sizeof(hit_charge) );
    memset( hit_ph, 0, sizeof(hit_ph) );
    memset( hit_startT, 0, sizeof(hit_startT) );
    memset( hit_endT, 0, sizeof(hit_endT) );
    memset( hit_rms, 0, sizeof(hit_rms) );
    memset( hit_electrons, 0, sizeof(hit_electrons) );

    Em_ch = 0;
    Em_e = 0;
    trk_e = 0;
    Emichel_e = 0;

    n_recoShowers = 0;
    memset( sh_direction_X, 0, sizeof(sh_direction_X) );
    memset( sh_direction_Y, 0, sizeof(sh_direction_Y) );
    memset( sh_direction_Z, 0, sizeof(sh_direction_Z) );
    memset( sh_start_X, 0, sizeof(sh_start_X) );
    memset( sh_start_Y, 0, sizeof(sh_start_Y) );
    memset( sh_start_Z, 0, sizeof(sh_start_Z) );
    memset( sh_energy, 0, sizeof(sh_energy) );
    memset( sh_MIPenergy, 0, sizeof(sh_MIPenergy) );
    memset( sh_dEdx, 0, sizeof(sh_dEdx) );
    memset( sh_bestplane, 0, sizeof(sh_bestplane) );
    memset( sh_length, 0, sizeof(sh_length) );

    n_flashes = 0;
    memset( flash_time, 0, sizeof(flash_time) );
    memset( flash_pe, 0, sizeof(flash_pe) );
    memset( flash_ycenter, 0, sizeof(flash_ycenter) );
    memset( flash_zcenter, 0, sizeof(flash_zcenter) );
    memset( flash_ywidth, 0, sizeof(flash_ywidth) );
    memset( flash_zwidth, 0, sizeof(flash_zwidth) );
    memset( flash_timewidth, 0, sizeof(flash_timewidth) );
    memset( flash_abstime, 0, sizeof(flash_abstime) );
    memset( flash_frame, 0, sizeof(flash_frame) );
    memset( flash_PE_ndk, 0, sizeof(flash_PE_ndk) );
    memset( flash_PE_Ar39, 0, sizeof(flash_PE_Ar39) );
  } // function NDKAna::reset

  //========================================================================
  DEFINE_ART_MODULE(NDKAna)

} // namespace DUNE
