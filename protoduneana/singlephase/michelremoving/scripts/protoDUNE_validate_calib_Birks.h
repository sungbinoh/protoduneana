//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep  9 15:17:32 2019 by ROOT version 6.16/00
// from TTree Event/Event Tree from Reco
// found on file: /dune/data2/users/apaudel/prod2_calib/r5387_sce.root
//////////////////////////////////////////////////////////

#ifndef protoDUNE_validate_calib_Birks_h
#define protoDUNE_validate_calib_Birks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.

class protoDUNE_validate_calib_Birks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Double_t        evttime;
   Int_t           run;
   Int_t           subrun;
   Int_t           year_month_date;
   Int_t           hour_min_sec;
   Int_t           cross_trks;
   Int_t           stopping_trks;
   Int_t           all_trks;
   Int_t           unbroken_trks;
   Float_t         trackthetaxz[30];   //[cross_trks]
   Float_t         trackthetayz[30];   //[cross_trks]
   Float_t         trkstartx[30];   //[cross_trks]
   Float_t         trkstarty[30];   //[cross_trks]
   Float_t         trkstartz[30];   //[cross_trks]
   Float_t         trkendx[30];   //[cross_trks]
   Float_t         trkendy[30];   //[cross_trks]
   Float_t         trkendz[30];   //[cross_trks]
   Float_t         trklen[30];   //[cross_trks]
   Float_t         peakT_max[30];   //[cross_trks]
   Float_t         peakT_min[30];   //[cross_trks]
   Int_t           TrkID[30];   //[cross_trks]
   Float_t         trkstartcosxyz[30][3];   //[cross_trks]
   Float_t         trkendcosxyz[30][3];   //[cross_trks]
   Int_t           ntrkhits[30][3];   //[cross_trks]
   Float_t         trkdqdx[30][3][3000];   //[cross_trks]
   Float_t         trkdedx[30][3][3000];   //[cross_trks]
   Float_t         trkresrange[30][3][3000];   //[cross_trks]
   Float_t         trkhitx[30][3][3000];   //[cross_trks]
   Float_t         trkhity[30][3][3000];   //[cross_trks]
   Float_t         trkhitz[30][3][3000];   //[cross_trks]
   Float_t         trkpitch[30][3][3000];   //[cross_trks]
   Float_t         dist_min[30];   //[cross_trks]
   Int_t           adjacent_hits[30];   //[cross_trks]
   Int_t           lastwire[30];   //[cross_trks]
   Float_t         lastpeakt[30];   //[cross_trks]
   Int_t           endtpc[30];   //[cross_trks]
   Float_t         true_trkstartx[30];   //[cross_trks]
   Float_t         true_trkstarty[30];   //[cross_trks]
   Float_t         true_trkstartz[30];   //[cross_trks]
   Float_t         true_trkendx[30];   //[cross_trks]
   Float_t         true_trkendy[30];   //[cross_trks]
   Float_t         true_trkendz[30];   //[cross_trks]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_run;   //!
   TBranch        *b_surbrun;   //!
   TBranch        *b_year_month_date;   //!
   TBranch        *b_hour_min_sec;   //!
   TBranch        *b_cross_trks;   //!
   TBranch        *b_stopping_trks;   //!
   TBranch        *b_all_trks;   //!
   TBranch        *b_unbroken_trks;   //!
   TBranch        *b_trackthetaxz;   //!
   TBranch        *b_trackthetayz;   //!
   TBranch        *b_trkstartx;   //!
   TBranch        *b_trkstarty;   //!
   TBranch        *b_trkstartz;   //!
   TBranch        *b_trkendx;   //!
   TBranch        *b_trkendy;   //!
   TBranch        *b_trkendz;   //!
   TBranch        *b_trklen;   //!
   TBranch        *b_peakT_max;   //!
   TBranch        *b_peakT_min;   //!
   TBranch        *b_TrkID;   //!
   TBranch        *b_trkstartcosxyz;   //!
   TBranch        *b_trkendcosxyz;   //!
   TBranch        *b_ntrkhits;   //!
   TBranch        *b_trkdqdx;   //!
   TBranch        *b_trkdedx;   //!
   TBranch        *b_trkresrange;   //!
   TBranch        *b_trkhitx;   //!
   TBranch        *b_trkhity;   //!
   TBranch        *b_trkhitz;   //!
   TBranch        *b_trkpitch;   //!
   TBranch        *b_dist_min;   //!
   TBranch        *b_adjacent_hits;   //!
   TBranch        *b_lastwire;   //!
   TBranch        *b_lastpeakt;   //!
   TBranch        *b_endtpc;   //!
   TBranch        *b_true_trkstartx;   //!
   TBranch        *b_true_trkstarty;   //!
   TBranch        *b_true_trkstartz;   //!
   TBranch        *b_true_trkendx;   //!
   TBranch        *b_true_trkendy;   //!
   TBranch        *b_true_trkendz;   //!

   protoDUNE_validate_calib_Birks(TTree *tree=0);
   virtual ~protoDUNE_validate_calib_Birks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int mn);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef protoDUNE_validate_calib_Birks_cxx
protoDUNE_validate_calib_Birks::protoDUNE_validate_calib_Birks(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/dune/data2/users/apaudel/prod2_calib/r5387_sce.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/dune/data2/users/apaudel/prod2_calib/r5387_sce.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/dune/data2/users/apaudel/prod2_calib/r5387_sce.root:/michelremoving");
      dir->GetObject("Event",tree);

   }
   Init(tree);
}

protoDUNE_validate_calib_Birks::~protoDUNE_validate_calib_Birks()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t protoDUNE_validate_calib_Birks::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t protoDUNE_validate_calib_Birks::LoadTree(Long64_t entry)
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

void protoDUNE_validate_calib_Birks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_surbrun);
   fChain->SetBranchAddress("year_month_date", &year_month_date, &b_year_month_date);
   fChain->SetBranchAddress("hour_min_sec", &hour_min_sec, &b_hour_min_sec);
   fChain->SetBranchAddress("cross_trks", &cross_trks, &b_cross_trks);
   fChain->SetBranchAddress("stopping_trks", &stopping_trks, &b_stopping_trks);
   fChain->SetBranchAddress("all_trks", &all_trks, &b_all_trks);
   fChain->SetBranchAddress("unbroken_trks", &unbroken_trks, &b_unbroken_trks);
   fChain->SetBranchAddress("trackthetaxz", trackthetaxz, &b_trackthetaxz);
   fChain->SetBranchAddress("trackthetayz", trackthetayz, &b_trackthetayz);
   fChain->SetBranchAddress("trkstartx", trkstartx, &b_trkstartx);
   fChain->SetBranchAddress("trkstarty", trkstarty, &b_trkstarty);
   fChain->SetBranchAddress("trkstartz", trkstartz, &b_trkstartz);
   fChain->SetBranchAddress("trkendx", trkendx, &b_trkendx);
   fChain->SetBranchAddress("trkendy", trkendy, &b_trkendy);
   fChain->SetBranchAddress("trkendz", trkendz, &b_trkendz);
   fChain->SetBranchAddress("trklen", trklen, &b_trklen);
   fChain->SetBranchAddress("peakT_max", peakT_max, &b_peakT_max);
   fChain->SetBranchAddress("peakT_min", peakT_min, &b_peakT_min);
   fChain->SetBranchAddress("TrkID", TrkID, &b_TrkID);
   fChain->SetBranchAddress("trkstartcosxyz", trkstartcosxyz, &b_trkstartcosxyz);
   fChain->SetBranchAddress("trkendcosxyz", trkendcosxyz, &b_trkendcosxyz);
   fChain->SetBranchAddress("ntrkhits", ntrkhits, &b_ntrkhits);
   fChain->SetBranchAddress("trkdqdx", trkdqdx, &b_trkdqdx);
   fChain->SetBranchAddress("trkdedx", trkdedx, &b_trkdedx);
   fChain->SetBranchAddress("trkresrange", trkresrange, &b_trkresrange);
   fChain->SetBranchAddress("trkhitx", trkhitx, &b_trkhitx);
   fChain->SetBranchAddress("trkhity", trkhity, &b_trkhity);
   fChain->SetBranchAddress("trkhitz", trkhitz, &b_trkhitz);
   fChain->SetBranchAddress("trkpitch", trkpitch, &b_trkpitch);
   fChain->SetBranchAddress("dist_min", dist_min, &b_dist_min);
   fChain->SetBranchAddress("adjacent_hits", adjacent_hits, &b_adjacent_hits);
   fChain->SetBranchAddress("lastwire", lastwire, &b_lastwire);
   fChain->SetBranchAddress("lastpeakt", lastpeakt, &b_lastpeakt);
   fChain->SetBranchAddress("endtpc", endtpc, &b_endtpc);
   TBranch *br = (TBranch*)fChain->GetListOfBranches()->FindObject("true_trkstartx");
   if (br){
     fChain->SetBranchAddress("true_trkstartx", true_trkstartx, &b_true_trkstartx);
     fChain->SetBranchAddress("true_trkstarty", true_trkstarty, &b_true_trkstarty);
     fChain->SetBranchAddress("true_trkstartz", true_trkstartz, &b_true_trkstartz);
     fChain->SetBranchAddress("true_trkendx", true_trkendx, &b_true_trkendx);
     fChain->SetBranchAddress("true_trkendy", true_trkendy, &b_true_trkendy);
     fChain->SetBranchAddress("true_trkendz", true_trkendz, &b_true_trkendz);
   }
   Notify();
}

Bool_t protoDUNE_validate_calib_Birks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void protoDUNE_validate_calib_Birks::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t protoDUNE_validate_calib_Birks::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef protoDUNE_validate_calib_Birks_cxx
