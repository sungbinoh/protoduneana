/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 23 17:37:10 2018 by ROOT version 6.13/01
// from TTree Event/Event Tree from Sim & Reco
// found on file: ../Data_reco_v2/ndk_test.root
//////////////////////////////////////////////////////////

#ifndef NDKAna_h
#define NDKAna_h

#include "TROOT.h"
//#include "/cvmfs/larsoft.opensciencegrid.org/products/root/v6_22_06a/Linux64bit+3.10-2.17-e19-p383b-debug/include/TROOT.h"

#include <TChain.h>
#include <TFile.h>
//#include "TH1F.h"
//#include "TH2F.h"
#include "TLorentzVector.h"
#include "TInterpreter.h"
#include <cstring>

#include "Limits.h"


// Header file for the classes stored in the TTree if any.
#include "vector"
using namespace std;
class NDKAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventNo;
   Int_t           runNo;
   Int_t           subRunNo;
   Double_t        MC_Ev;
   Int_t           MC_cc;
   Double_t        MC_Q2;
   Int_t           MC_nuPDG;
   Double_t        MC_hit_nucleon;
   Int_t           mcgenie_npart;
   Int_t           mcgenie_id[MAX_GEN_PARTS];   //[mcgenie_npart]
   Int_t           mcgenie_fate[MAX_GEN_PARTS];   //[mcgenie_npart]
   Int_t           mcgenie_statusCode[MAX_GEN_PARTS];   //[mcgenie_npart]
   Int_t           mcgenie_pdg[MAX_GEN_PARTS];   //[mcgenie_npart]
   Int_t           mcgenie_mother[MAX_GEN_PARTS];   //[mcgenie_npart]
   Double_t        mcgenie_startMomentum[MAX_GEN_PARTS][4];   //[mcgenie_npart]
   Double_t        mcgenie_endMomentum[MAX_GEN_PARTS][4];   //[mcgenie_npart]
   Double_t        mc_vertex[4];
   Int_t           mc_npart;
   Int_t           mc_id[MAX_MC_PARTS];   //[mc_npart]
   Int_t           mc_pdg[MAX_MC_PARTS];   //[mc_npart]
   Int_t           mc_statusCode[MAX_MC_PARTS];   //[mc_npart]
   Int_t           mc_mother[MAX_MC_PARTS];   //[mc_npart]
   Double_t        mc_startXYZT[MAX_MC_PARTS][4];   //[mc_npart]
   Double_t        mc_endXYZT[MAX_MC_PARTS][4];   //[mc_npart]
   int             mc_npoints[MAX_MC_PARTS];
   double          mc_xyz[MAX_MC_PARTS][MAX_MC_TRAJ_PTS][3];
   Double_t        mc_startMomentum[MAX_MC_PARTS][4];   //[mc_npart]
   Double_t        mc_endMomentum[MAX_MC_PARTS][4];   //[mc_npart]
   Double_t        mc_Prange[MAX_MC_PARTS];   //[mc_npart]
   Double_t        mc_truthlength[MAX_MC_PARTS];   //[mc_npart]
   Double_t        MC_length_Kaons[MAX_MC_PARTS];
   Double_t        MC_length_Muons[MAX_MC_PARTS];
   Double_t        MC_length_Pions[MAX_MC_PARTS];
   vector<string>  *mc_process;
   vector<string>  *mc_Endprocess;

   Int_t           n_vertices;
   Double_t        vertex[MAX_RECO_VERTICES][4];   //[n_vertices]
   Int_t           n_reco_tracks;
   Int_t           n_decayVtx;
   Double_t        decayVtx[MAX_RECO_VERTICES][3];   //[n_decayVtx]
   Double_t        track_vtx[MAX_RECO_TRACKS][4];   //[n_reco_tracks]
   Double_t        track_vtxDir[MAX_RECO_TRACKS][3];   //[n_reco_tracks]
   Double_t        track_end[MAX_RECO_TRACKS][4];   //[n_reco_tracks]
   Int_t           track_isContained[MAX_RECO_TRACKS];   //[n_reco_tracks]
   Double_t        track_length[MAX_RECO_TRACKS];   //[n_reco_tracks]
   Double_t        track_length_Kaons[MAX_RECO_TRACKS];
   Double_t        track_length_Muons[MAX_RECO_TRACKS];
   Double_t        track_length_Pions[MAX_RECO_TRACKS];
   Double_t        track_PIDA[MAX_RECO_TRACKS][3];   //[n_reco_tracks]
   Int_t           track_PID_pdg[MAX_RECO_TRACKS][3];   //[n_reco_tracks]
   Double_t        track_KE[MAX_RECO_TRACKS][3];   //[n_reco_tracks]
   Double_t        track_Prange[MAX_RECO_TRACKS];   //[n_reco_tracks]
   Int_t           track_bestplane[MAX_RECO_TRACKS];   //[n_reco_tracks]
   Int_t           n_cal_points[MAX_RECO_TRACKS];
   Double_t        track_dQ_dx[MAX_RECO_TRACKS][MAX_CALO_PTS];   //[n_cal_points]
   Double_t        track_dE_dx[MAX_RECO_TRACKS][MAX_CALO_PTS];   //[n_cal_points]
   Double_t        track_range[MAX_RECO_TRACKS][MAX_CALO_PTS];   //[n_cal_points]
   Double_t        track_pitch[MAX_RECO_TRACKS][MAX_CALO_PTS];   //[n_cal_points]
   Int_t           n_cal_points_byplane[MAX_RECO_TRACKS][3];
   Double_t        track_dQ_dx_byplane[MAX_RECO_TRACKS][3][MAX_CALO_PTS];   //[n_cal_points]
   Double_t        track_dE_dx_byplane[MAX_RECO_TRACKS][3][MAX_CALO_PTS];   //[n_cal_points]
   Double_t        track_range_byplane[MAX_RECO_TRACKS][3][MAX_CALO_PTS];   //[n_cal_points]
   Double_t        track_pitch_byplane[MAX_RECO_TRACKS][3][MAX_CALO_PTS];   //[n_cal_points]
   Double_t        track_calo_xyz_byplane[MAX_RECO_TRACKS][3][MAX_CALO_PTS][3];   //[n_reco_tracks]
   Double_t        track_complet[MAX_RECO_TRACKS];   //[n_reco_tracks]
   Double_t        track_Efrac[MAX_RECO_TRACKS];   //[n_reco_tracks]
   Int_t           track_mcID[MAX_RECO_TRACKS];   //[n_reco_tracks]
   Int_t           track_mcPDG[MAX_RECO_TRACKS];   //[n_reco_tracks]
    int    n_track_points[MAX_RECO_TRACKS];
    double track_point_xyz[MAX_RECO_TRACKS][MAX_CALO_PTS][3];
   Double_t        Em_ch;
   Double_t        Em_e;
   Double_t        trk_e;
   Double_t        Emichel_e;
   Int_t           n_showers;
   Double_t        sh_direction_X[MAX_SHOWERS];   //[n_showers]
   Double_t        sh_direction_Y[MAX_SHOWERS];   //[n_showers]
   Double_t        sh_direction_Z[MAX_SHOWERS];   //[n_showers]
   Double_t        sh_start_X[MAX_SHOWERS];   //[n_showers]
   Double_t        sh_start_Y[MAX_SHOWERS];   //[n_showers]
   Double_t        sh_start_Z[MAX_SHOWERS];   //[n_showers]
   Double_t        sh_energy[MAX_SHOWERS][3];   //[n_showers]
   Double_t        sh_MIPenergy[MAX_SHOWERS][3];   //[n_showers]
   Double_t        sh_dEdx[MAX_SHOWERS][3];   //[n_showers]
   Int_t           sh_bestplane[MAX_SHOWERS];   //[n_showers]
   Double_t        sh_length[MAX_SHOWERS];   //[n_showers]
   Int_t           n_flashes;
   Double_t        flash_time[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_pe[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_ycenter[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_zcenter[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_ywidth[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_zwidth[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_timewidth[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_abstime[MAX_FLASHES];   //[n_flashes]
   Int_t           flash_frame[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_PE_ndk[MAX_FLASHES];   //[n_flashes]
   Double_t        flash_PE_Ar39[MAX_FLASHES];   //[n_flashes]

   // List of branches
   TBranch        *b_eventNo;   //!
   TBranch        *b_runNo;   //!
   TBranch        *b_subRunNo;   //!
   TBranch        *b_MC_Ev;   //!
   TBranch        *b_MC_cc;   //!
   TBranch        *b_MC_Q2;   //!
   TBranch        *b_MC_nuPDG;   //!
   TBranch        *b_MC_hit_nucleon;   //!
   TBranch        *b_mcgenie_npart;   //!
   TBranch        *b_mcgenie_id;   //!
   TBranch        *b_mcgenie_fate;   //!
   TBranch        *b_mcgenie_statusCode;   //!
   TBranch        *b_mcgenie_pdg;   //!
   TBranch        *b_mcgenie_mother;   //!
   TBranch        *b_mcgenie_startMomentum;   //!
   TBranch        *b_mcgenie_endMomentum;   //!
   TBranch        *b_mc_vertex;   //!
   TBranch        *b_mc_npart;   //!
   TBranch        *b_mc_id;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_statusCode;   //!
   TBranch        *b_mc_mother;   //!
   TBranch        *b_mc_startXYZT;   //!
   TBranch        *b_mc_endXYZT;   //!
   TBranch        *b_mc_npoints;   //!
   TBranch        *b_mc_xyz;   //!
   TBranch        *b_mc_startMomentum;   //!
   TBranch        *b_mc_endMomentum;   //!
   TBranch        *b_mc_Prange;   //!
   TBranch        *b_mc_truthlength;   //!
   TBranch        *b_MC_length_Kaons; //!
   TBranch        *b_MC_length_Muons; //!
   TBranch        *b_MC_length_Pions; //!
  // TBranch        *b_mc_process;   //!
  // TBranch        *b_mc_Endprocess;   //!
   TBranch        *b_n_vertices;   //!
   TBranch        *b_vertex;   //!
   TBranch        *b_n_reco_tracks;   //!
   TBranch        *b_n_decayVtx;   //!
   TBranch        *b_decayVtx;   //!
   TBranch        *b_track_vtx;   //!
   TBranch        *b_track_vtxDir;   //!
   TBranch        *b_track_end;   //!
   TBranch        *b_track_isContained;   //!
   TBranch        *b_track_length;   //!
   TBranch        *b_track_length_Kaons; //!
   TBranch        *b_track_length_Muons; //!
   TBranch        *b_track_length_Pions; //!
   TBranch        *b_track_PIDA;   //!
   TBranch        *b_track_PID_pdg;   //!
   TBranch        *b_track_KE;   //!
   TBranch        *b_track_Prange;   //!
   TBranch        *b_track_bestplane;   //!
   TBranch        *b_n_cal_points;   //!
   TBranch        *b_track_dQ_dx;   //!
   TBranch        *b_track_dE_dx;   //!
   TBranch        *b_track_range;   //!
   TBranch        *b_track_pitch;   //!
   TBranch        *b_n_cal_points_byplane;   //!
   TBranch        *b_track_dQ_dx_byplane;
   TBranch        *b_track_dE_dx_byplane;
   TBranch        *b_track_range_byplane;
   TBranch        *b_track_pitch_byplane;
   TBranch        *b_track_calo_xyz_byplane;   //!
   TBranch        *b_track_complet;   //!
   TBranch        *b_track_Efrac;   //!
   TBranch        *b_track_mcID;   //!
   TBranch        *b_track_mcPDG;   //!
  // TBranch        *b_n_track_points;   //!
//   TBranch        *b_track_point_xyz;   //!
   TBranch        *b_Em_ch;   //!
   TBranch        *b_Em_e;   //!
   TBranch        *b_trk_e;   //!
   TBranch        *b_Emichel_e;   //!
   TBranch        *b_n_showers;   //!
   TBranch        *b_sh_direction_X;   //!
   TBranch        *b_sh_direction_Y;   //!
   TBranch        *b_sh_direction_Z;   //!
   TBranch        *b_sh_start_X;   //!
   TBranch        *b_sh_start_Y;   //!
   TBranch        *b_sh_start_Z;   //!
   TBranch        *b_sh_energy;   //!
   TBranch        *b_sh_MIPenergy;   //!
   TBranch        *b_sh_dEdx;   //!
   TBranch        *b_sh_bestplane;   //!
   TBranch        *b_sh_length;   //!
   TBranch        *b_n_flashes;   //!
   TBranch        *b_flash_time;   //!
   TBranch        *b_flash_pe;   //!
   TBranch        *b_flash_ycenter;   //!
   TBranch        *b_flash_zcenter;   //!
   TBranch        *b_flash_ywidth;   //!
   TBranch        *b_flash_zwidth;   //!
   TBranch        *b_flash_timewidth;   //!
   TBranch        *b_flash_abstime;   //!
   TBranch        *b_flash_frame;   //!
   TBranch        *b_flash_PE_ndk;   //!
   TBranch        *b_flash_PE_Ar39;   //!

   NDKAna(TTree *tree=0);
   NDKAna(const char *fnames);
   virtual ~NDKAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void     Reset();
};

#endif

 

