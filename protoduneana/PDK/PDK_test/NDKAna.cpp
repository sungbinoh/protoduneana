#define NDKAna_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "NDKAna.h"
#include <TRegexp.h>




NDKAna::NDKAna(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Data_reco_v2/ndk_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../Data_reco_v2/ndk_test.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../Data_reco_v2/ndk_test.root:/NDKAna");
      dir->GetObject("Event",tree);

   }
   Init(tree);
}

NDKAna::NDKAna(const char *fnames)
{
    TString names(fnames);
    TRegexp pat("[^ ]*\\.root");

    TChain* tree = new TChain("NDKAna/Event");

    TSubString subs(names(pat));
    while ( !subs.IsNull() ) {
        tree->Add(TString(subs));
        names.Remove(subs.Start(), subs.Length());
        subs = names(pat);
    }

    Init(tree);
}

NDKAna::~NDKAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NDKAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NDKAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NDKAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mc_process = 0;
   mc_Endprocess = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNo", &eventNo, &b_eventNo);
   fChain->SetBranchAddress("runNo", &runNo, &b_runNo);
   fChain->SetBranchAddress("subRunNo", &subRunNo, &b_subRunNo);
   fChain->SetBranchAddress("MC_Ev", &MC_Ev, &b_MC_Ev);
   fChain->SetBranchAddress("MC_cc", &MC_cc, &b_MC_cc);
   fChain->SetBranchAddress("MC_Q2", &MC_Q2, &b_MC_Q2);
   fChain->SetBranchAddress("MC_nuPDG", &MC_nuPDG, &b_MC_nuPDG);
   fChain->SetBranchAddress("MC_hit_nucleon", &MC_hit_nucleon, &b_MC_hit_nucleon);
   fChain->SetBranchAddress("mcgenie_npart", &mcgenie_npart, &b_mcgenie_npart);
   fChain->SetBranchAddress("mcgenie_id", mcgenie_id, &b_mcgenie_id);
   fChain->SetBranchAddress("mcgenie_fate", mcgenie_fate, &b_mcgenie_fate);
   fChain->SetBranchAddress("mcgenie_statusCode", mcgenie_statusCode, &b_mcgenie_statusCode);
   fChain->SetBranchAddress("mcgenie_pdg", mcgenie_pdg, &b_mcgenie_pdg);
   fChain->SetBranchAddress("mcgenie_mother", mcgenie_mother, &b_mcgenie_mother);
   fChain->SetBranchAddress("mcgenie_startMomentum", mcgenie_startMomentum, &b_mcgenie_startMomentum);
   fChain->SetBranchAddress("mcgenie_endMomentum", mcgenie_endMomentum, &b_mcgenie_endMomentum);
   fChain->SetBranchAddress("mc_vertex", mc_vertex, &b_mc_vertex);
   fChain->SetBranchAddress("mc_npart", &mc_npart, &b_mc_npart);
   fChain->SetBranchAddress("mc_id", mc_id, &b_mc_id);
   fChain->SetBranchAddress("mc_pdg", mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_statusCode", mc_statusCode, &b_mc_statusCode);
   fChain->SetBranchAddress("mc_mother", mc_mother, &b_mc_mother);
   fChain->SetBranchAddress("mc_startXYZT", mc_startXYZT, &b_mc_startXYZT);
   fChain->SetBranchAddress("mc_endXYZT", mc_endXYZT, &b_mc_endXYZT);
   fChain->SetBranchAddress("mc_npoints", mc_npoints, &b_mc_npoints);
   fChain->SetBranchAddress("mc_xyz", mc_xyz, &b_mc_xyz);
   fChain->SetBranchAddress("mc_startMomentum", mc_startMomentum, &b_mc_startMomentum);
   fChain->SetBranchAddress("mc_endMomentum", mc_endMomentum, &b_mc_endMomentum);
   fChain->SetBranchAddress("mc_Prange", mc_Prange, &b_mc_Prange);
   fChain->SetBranchAddress("mc_truthlength", mc_truthlength, &b_mc_truthlength);
   fChain->SetBranchAddress("MC_length_Kaons", MC_length_Kaons, &b_MC_length_Kaons);
   fChain->SetBranchAddress("MC_length_Muons", MC_length_Muons, &b_MC_length_Muons);
   fChain->SetBranchAddress("MC_length_Pions", MC_length_Pions, &b_MC_length_Pions);



 //  fChain->SetBranchAddress("mc_process", &mc_process, &b_mc_process);
  // fChain->SetBranchAddress("mc_Endprocess", &mc_Endprocess, &b_mc_Endprocess);
   fChain->SetBranchAddress("n_vertices", &n_vertices, &b_n_vertices);
   fChain->SetBranchAddress("vertex", vertex, &b_vertex);
   fChain->SetBranchAddress("n_reco_tracks", &n_reco_tracks, &b_n_reco_tracks);
   fChain->SetBranchAddress("n_decayVtx", &n_decayVtx, &b_n_decayVtx);
   fChain->SetBranchAddress("decayVtx", &decayVtx, &b_decayVtx);
   fChain->SetBranchAddress("track_vtx", track_vtx, &b_track_vtx);
   fChain->SetBranchAddress("track_vtxDir", track_vtxDir, &b_track_vtxDir);
   fChain->SetBranchAddress("track_end", track_end, &b_track_end);
   fChain->SetBranchAddress("track_isContained", track_isContained, &b_track_isContained);
   fChain->SetBranchAddress("track_length", track_length, &b_track_length);
   fChain->SetBranchAddress("track_length_Kaons", track_length_Kaons, &b_track_length_Kaons);
   fChain->SetBranchAddress("track_length_Muons", track_length_Muons, &b_track_length_Muons);
   fChain->SetBranchAddress("track_length_Pions", track_length_Pions, &b_track_length_Pions);



   fChain->SetBranchAddress("track_PIDA", track_PIDA, &b_track_PIDA);
   fChain->SetBranchAddress("track_PID_pdg", track_PID_pdg, &b_track_PID_pdg);
   fChain->SetBranchAddress("track_KE", track_KE, &b_track_KE);
   fChain->SetBranchAddress("track_Prange", track_Prange, &b_track_Prange);
   fChain->SetBranchAddress("track_bestplane", track_bestplane, &b_track_bestplane);
   fChain->SetBranchAddress("n_cal_points", n_cal_points, &b_n_cal_points);
   fChain->SetBranchAddress("track_dQ_dx", track_dQ_dx, &b_track_dQ_dx);
   fChain->SetBranchAddress("track_dE_dx", track_dE_dx, &b_track_dE_dx);
   fChain->SetBranchAddress("track_range", track_range, &b_track_range);
   fChain->SetBranchAddress("track_pitch", track_pitch, &b_track_pitch);
   fChain->SetBranchAddress("n_cal_points_byplane", n_cal_points_byplane, &b_n_cal_points_byplane);
   fChain->SetBranchAddress("track_dQ_dx_byplane", track_dQ_dx_byplane, &b_track_dQ_dx_byplane);
   fChain->SetBranchAddress("track_dE_dx_byplane", track_dE_dx_byplane, &b_track_dE_dx_byplane);
   fChain->SetBranchAddress("track_range_byplane", track_range_byplane, &b_track_range_byplane);
   fChain->SetBranchAddress("track_pitch_byplane", track_pitch_byplane, &b_track_pitch_byplane);
   fChain->SetBranchAddress("track_calo_xyz_byplane", track_calo_xyz_byplane, &b_track_calo_xyz_byplane);
   fChain->SetBranchAddress("track_complet", track_complet, &b_track_complet);
   fChain->SetBranchAddress("track_Efrac", track_Efrac, &b_track_Efrac);
   fChain->SetBranchAddress("track_mcID", track_mcID, &b_track_mcID);
   fChain->SetBranchAddress("track_mcPDG", track_mcPDG, &b_track_mcPDG);
  // fChain->SetBranchAddress("n_track_points", n_track_points, &b_n_track_points);
  // fChain->SetBranchAddress("track_point_xyz", track_point_xyz, &b_track_point_xyz);
   fChain->SetBranchAddress("Em_ch", &Em_ch, &b_Em_ch);
   fChain->SetBranchAddress("Em_e", &Em_e, &b_Em_e);
   fChain->SetBranchAddress("trk_e", &trk_e, &b_trk_e);
   fChain->SetBranchAddress("Emichel_e", &Emichel_e, &b_Emichel_e);
   fChain->SetBranchAddress("n_showers", &n_showers, &b_n_showers);
   fChain->SetBranchAddress("sh_direction_X", sh_direction_X, &b_sh_direction_X);
   fChain->SetBranchAddress("sh_direction_Y", sh_direction_Y, &b_sh_direction_Y);
   fChain->SetBranchAddress("sh_direction_Z", sh_direction_Z, &b_sh_direction_Z);
   fChain->SetBranchAddress("sh_start_X", sh_start_X, &b_sh_start_X);
   fChain->SetBranchAddress("sh_start_Y", sh_start_Y, &b_sh_start_Y);
   fChain->SetBranchAddress("sh_start_Z", sh_start_Z, &b_sh_start_Z);
   fChain->SetBranchAddress("sh_energy", sh_energy, &b_sh_energy);
   fChain->SetBranchAddress("sh_MIPenergy", sh_MIPenergy, &b_sh_MIPenergy);
   fChain->SetBranchAddress("sh_dEdx", sh_dEdx, &b_sh_dEdx);
   fChain->SetBranchAddress("sh_bestplane", sh_bestplane, &b_sh_bestplane);
   fChain->SetBranchAddress("sh_length", sh_length, &b_sh_length);
   fChain->SetBranchAddress("n_flashes", &n_flashes, &b_n_flashes);
   fChain->SetBranchAddress("flash_time", &flash_time, &b_flash_time);
   fChain->SetBranchAddress("flash_pe", &flash_pe, &b_flash_pe);
   fChain->SetBranchAddress("flash_ycenter", &flash_ycenter, &b_flash_ycenter);
   fChain->SetBranchAddress("flash_zcenter", &flash_zcenter, &b_flash_zcenter);
   fChain->SetBranchAddress("flash_ywidth", &flash_ywidth, &b_flash_ywidth);
   fChain->SetBranchAddress("flash_zwidth", &flash_zwidth, &b_flash_zwidth);
   fChain->SetBranchAddress("flash_timewidth", &flash_timewidth, &b_flash_timewidth);
   fChain->SetBranchAddress("flash_abstime", &flash_abstime, &b_flash_abstime);
   fChain->SetBranchAddress("flash_frame", &flash_frame, &b_flash_frame);
   fChain->SetBranchAddress("flash_PE_ndk", &flash_PE_ndk, &b_flash_PE_ndk);
   fChain->SetBranchAddress("flash_PE_Ar39", &flash_PE_Ar39, &b_flash_PE_Ar39);
   Notify();
}

Bool_t NDKAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NDKAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NDKAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#define RESET(i) std::memset(i, 0, sizeof(i))
#define RESETSINGLE(i) i = 0
                                                                                                     

void NDKAna::Reset()
{
    size_t howmuch = (const char*)flash_PE_Ar39 - (const char*)&eventNo
        + sizeof(flash_PE_Ar39);
    std::memset(&eventNo, 0, howmuch);
    /*
    RESETSINGLE(eventNo);
    RESETSINGLE(runNo);
    RESETSINGLE(subRunNo);
    RESETSINGLE(MC_Ev);
    RESETSINGLE(MC_cc);
    RESETSINGLE(MC_Q2);
    RESETSINGLE(MC_nuPDG);
    RESETSINGLE(MC_hit_nucleon);
    RESETSINGLE(mcgenie_npart);
    RESET(mcgenie_id);   //[mcgenie_npart]
    RESET(mcgenie_fate);   //[mcgenie_npart]
    RESET(mcgenie_statusCode);   //[mcgenie_npart]
    RESET(mcgenie_pdg);   //[mcgenie_npart]
    RESET(mcgenie_mother);   //[mcgenie_npart]
    RESET(mcgenie_startMomentum);   //[mcgenie_npart]
    RESET(mcgenie_endMomentum);   //[mcgenie_npart]
    RESET(mc_vertex);
    RESETSINGLE(mc_npart);
    RESET(mc_id);
    RESET(mc_pdg);
    RESET(mc_statusCode);
    RESET(mc_mother);
    RESET(mc_startXYZT);
    RESET(mc_endXYZT);
    RESET(mc_npoints);
    RESET(mc_xyz);
    RESET(mc_startMomentum);
    RESET(mc_endMomentum);
    RESET(mc_Prange);
    RESET(mc_truthlength);

    RESETSINGLE(n_vertices);
    RESET(vertex);
    RESETSINGLE(n_reco_tracks);
    RESETSINGLE(n_decayVtx);
    RESET(decayVtx);
    RESET(track_vtx);
    RESET(track_vtxDir);
    RESET(track_end);
    RESET(track_isContained);
    RESET(track_length);
    RESET(track_PIDA);
    RESET(track_PID_pdg);
    RESET(track_KE);
    RESET(track_Prange);
    RESET(track_bestplane);
    RESET(n_cal_points);
    RESET(track_dQ_dx);
    RESET(track_dE_dx);
    RESET(track_range);
    RESET(track_pitch);
    RESET(n_cal_points_byplane);
    RESET(track_dQ_dx_byplane);
    RESET(track_dE_dx_byplane);
    RESET(track_range_byplane);
    RESET(track_pitch_byplane);
    RESET(track_calo_xyz_byplane);
    RESET(track_complet);
    RESET(track_Efrac);
    RESET(track_mcID);
    RESET(track_mcPDG);
    RESETSINGLE(Em_ch);
    RESETSINGLE(Em_e);
    RESETSINGLE(trk_e);
    RESETSINGLE(Emichel_e);
    RESETSINGLE(n_showers);
    RESET(sh_direction_X);
    RESET(sh_direction_Y);
    RESET(sh_direction_Z);
    RESET(sh_start_X);
    RESET(sh_start_Y);
    RESET(sh_start_Z);
    RESET(sh_energy);
    RESET(sh_MIPenergy);
    RESET(sh_dEdx);
    RESET(sh_bestplane);
    RESET(sh_length);
    RESETSINGLE(n_flashes);
    RESET(flash_time);
    RESET(flash_pe);
    RESET(flash_ycenter);
    RESET(flash_zcenter);
    RESET(flash_ywidth);
    RESET(flash_zwidth);
    RESET(flash_timewidth);
    RESET(flash_abstime);
    RESET(flash_frame);
    RESET(flash_PE_ndk);
    RESET(flash_PE_Ar39);
*/
}


#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void NDKAna::Loop()
{
//   In a ROOT session, you can do:
//      root> .L NDKAna.C
//      root> NDKAna t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
