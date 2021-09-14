//===================================================
// Template for log-likelihood PID
// Based on C. Marshall's idea
// see: https://indico.fnal.gov/event/18097/contribution/0/material/slides/0.pdf
// Skol...
//
// Modified to deal with reconstructed distance rather than calo-points/wire number
//===================================================
#include <iostream>
#include <algorithm>
#include "NDKAna.h"
#include "include.h"
using namespace std;


int Template_PIDloglike_range(  const char *filename, const char *outfile)
{

    //==================================================
    // Cut stats
    unsigned int cuts_stats[10] = {};

    enum cut_t {
	kNoCut,
	kTwoTracks,
	kShortCaloPts,
	kCommonVtx,
	kNoPi0,
	kNoLongProton,
	kShortNoLepton,
	kShortTooShort
    };

    const char* cut_names[] = {
	"No Cut POOP",
	"Two Tracks",
	"Calo pts in short",
	"Common Vertex",
	"No pi0",
	"Long is not P",
	"Short is no lepton/pion",
	"Too Short"
    };


    //TFile *f_data = new TFile(filename,"READ");

  TFile *f = new TFile(outfile,"RECREATE");
  // TH2D *h_dEdx = new TH2D("h_dEdx",";Distance from end vertex [cm];dE/dx",50,0,50,50,0,50);
  // TH2D *h_dEdx_rev = new TH2D("h_dEdx_rev",";Distance from opposite end [cm];dE/dx",50,0,50,50,0,50);
  TH2D *h_dEdx = new TH2D("h_dEdx",";Distance from end vertex [cm];dE/dx",50,0,30,50,0,50);
  TH2D *h_dEdx_rev = new TH2D("h_dEdx_rev",";Distance from opposite end [cm];dE/dx",50,0,30,50,0,50);
  
  TH1F *h_track_length_Kaons = new TH1F("h_track_length_Kaons", "Reco Track Length Kaons", 100, 0, 80);
  TH1F *h_track_length_Muons = new TH1F("h_track_length_Muons", "Reco Track Length Muons", 100, 0, 80); 
  TH1F *h_track_length_Pions = new TH1F("h_track_length_Pions", "Reco Track Length Pions", 100, 0, 80); 
  TH1F *h_track_length = new TH1F("h_track_length", "Reco Track Length", 100, 0, 80);
  TH1F *h_MC_length_Kaons = new TH1F("h_MC_length_Kaons","True Track Length Kaons", 100, 0, 80);
  TH1F *h_MC_length_Pions = new TH1F("h_MC_length_Pions","True Track Length Pions", 100, 0, 80);
  TH1F *h_MC_length_Muons = new TH1F("h_MC_length_Muons","True Track Length Muons", 100, 0, 80);
  TH1F *h_MC_truthlength = new TH1F("h_MC_truthlength","True Track Lengths", 100, 0, 80);



  TH1F *h_mc_pdg = new TH1F("h_mc_pdg", "True PDG", 100, 0, 2300);
  TH1F *h_track_PID_pdg = new TH1F("h_track_PID_pdg", "Reco PDG", 100, 0, 2300);

//========================================================
  TChain *tree = new TChain("NDKAna/Event");
  tree->Add(filename);
  NDKAna *signal = new NDKAna(tree);
  int n_entries = tree->GetEntries();

  // disable unnecessary branches
  vector<string> allowed ({"n_reco_tracks"
	      , "track_length"
	      , "n_cal_points"
	      , "track_vtx"
	      , "track_end"
	      , "mc_npart"
	      , "mc_mother"
	      , "mc_pdg"
	      , "track_mcPDG"
	      , "track_PID_pdg"
              , "track_range"
	      , "track_dE_dx"
              , "MC_length_Kaons"
              , "MC_length_Pions"
              , "MC_length_Muons"
              , "mc_truthlength"
              , "track_length_Muons"
              , "track_length_Kaons"
              , "track_length_Pions" });

  tree->SetBranchStatus("*", 0);
  for (auto it = allowed.begin(); it != allowed.end(); it++)
      tree->SetBranchStatus(it->c_str(), 1);

  // loop over the entries
  cout<<"..this many entries "<<n_entries<<endl;

//  h_track_length_Kaons -> Fill( track_length_Kaons );
  

  int fiftieth = n_entries / 50;
  cout<<"|                                                  |\r|";
  for( int i=0; i<n_entries; ++i){

    if ((i+1)%fiftieth == 0) {
      cout<<"-";
      cout.flush();
    }

     signal->GetEntry(i);
//h_track_length_Kaons -> Fill( track_length_Kaons );
 
    cuts_stats[kNoCut]++;

     //==========================================================================================
     // RECO INFO
     //==========================================================================================
//h_track_length_Kaons -> Fill(signal->track_length);
     // look at events with 2 tracks only!
//     if( signal->n_reco_tracks != 2 ) continue;
     cuts_stats[kTwoTracks]++;

     //================================================
     // longest and shortest tracks
     //================================================
     float length_tmp = -999.0;
     int longest_trk_idx =-1;
     float short_tmp = 999999;
     int shortest_trk_idx =-1;
     //int longest_vtx_ID =-999;
     for( int ii=0; ii<signal->n_reco_tracks; ++ii){
        h_track_length_Kaons -> Fill(signal->track_length_Kaons[ii]);
        h_track_length_Muons -> Fill(signal->track_length_Muons[ii]);
        h_track_length_Pions -> Fill(signal->track_length_Pions[ii]);
        h_track_length -> Fill(signal->track_length[ii]);
        
for (int j = 0; j<1; ++j){
        h_track_PID_pdg -> Fill(signal->track_PID_pdg[ii][j]);
}
	if( signal->track_length[ii] > length_tmp ){
          length_tmp = signal->track_length[ii];
          longest_trk_idx = ii;
        }
        if( signal->track_length[ii] < short_tmp){
           short_tmp = signal->track_length[ii];
           shortest_trk_idx = ii;
        }
     }
     if( shortest_trk_idx  == longest_trk_idx ) continue; //this should never happened though

     if ( signal -> n_cal_points[shortest_trk_idx] == 0 ) continue;
     cuts_stats[kShortCaloPts]++;

     // get tracks' vertices and ends
     TVector3 sh_trk_vtx(signal->track_vtx[shortest_trk_idx][0],
			 signal->track_vtx[shortest_trk_idx][1],
			 signal->track_vtx[shortest_trk_idx][2]);
     TVector3 sh_trk_end(signal->track_end[shortest_trk_idx][0],
			 signal->track_end[shortest_trk_idx][1],
			 signal->track_end[shortest_trk_idx][2]);
     TVector3 lo_trk_vtx(signal->track_vtx[longest_trk_idx][0],
			 signal->track_vtx[longest_trk_idx][1],
			 signal->track_vtx[longest_trk_idx][2]);
     TVector3 lo_trk_end(signal->track_end[longest_trk_idx][0],
			 signal->track_end[longest_trk_idx][1],
			 signal->track_end[longest_trk_idx][2]);

     double min_start = std::min((sh_trk_vtx-lo_trk_vtx).Mag(), (sh_trk_end-lo_trk_vtx).Mag() );
     double min_end = std::min((sh_trk_vtx-lo_trk_end).Mag(), (sh_trk_end-lo_trk_end).Mag() );

     //==== Cut on common vertex ====
     // the 2 tracks should have a common vertex! cut on 5 cm distance.
     if (std::min(min_start, min_end) > 5.) continue;
     cuts_stats[kCommonVtx]++;


     //========================================
     // Preselection
     //========================================

     //==== Cut on true primary pi0 - mother id == 0; affects only background ====
     int has_pi0 = 0;
     for (int j = 0; j < signal->mc_npart; j++) {

       h_mc_pdg -> Fill(signal->mc_pdg[j]);
       h_MC_length_Kaons -> Fill(signal->MC_length_Kaons[j]);
       h_MC_length_Pions -> Fill(signal->MC_length_Pions[j]);
       h_MC_length_Muons -> Fill(signal->MC_length_Muons[j]);
       h_MC_truthlength -> Fill(signal->mc_truthlength[j]);




	 if (signal->mc_mother[j] > 0) break; // not a primary particle
	 if (signal->mc_pdg[j] == 111) {// this is pi0
	     has_pi0 = 1;
	     break;
	 }
     }
     if (has_pi0) continue;
     cuts_stats[kNoPi0]++;

     //==== cut on longer track true ID - assuming we can tell it from proton ====
     if( abs(signal->track_mcPDG[longest_trk_idx]) >= 2212 )
	 //== 2212 ) // || signal->track_mcPDG[shortest_trk_idx] < 321 )
	 continue;
     cuts_stats[kNoLongProton]++;

     //==== cut on shorter track - either limit it to K and heavier, or allow muons, and pis
     if ( abs (signal->track_mcPDG[shortest_trk_idx]) < 321 ) continue;
     cuts_stats[kShortNoLepton]++;

     //========================================
     // determine short track direction
     //========================================
     int sh_trk_dir =0;
     int lo_trk_dir =0;
     //check muon and short track directions & common vertex?
     // dir = 1
     // long track pointing out of the vertex, short track pointing
     // towards the common vertex
     // dir =-1
     // vice versa
     if( min_start < min_end ){
       lo_trk_dir =1;
       std::cout<<"lo_trk_dir" << lo_trk_dir << endl;
       
	if( (sh_trk_vtx-lo_trk_vtx).Mag() < (sh_trk_end - lo_trk_vtx).Mag() ) sh_trk_dir = -1;
       else sh_trk_dir = 1;
     }
     else{
       lo_trk_dir =-1;
       std::cout<<"lo_trk_dir" << lo_trk_dir << endl;
      
 if( (sh_trk_vtx-lo_trk_end).Mag() < (sh_trk_end - lo_trk_end).Mag() ) sh_trk_dir = -1;
       else sh_trk_dir = 1;
       std::cout<<"lo_trk_dir" << lo_trk_dir << endl;

     }

     // 1 if 1st calo point near vertex and last near the end of the track, 0 otherwise
     int npoints = signal->n_cal_points[shortest_trk_idx];
     // 1st calo point always at the end of reco track!
     // Residual range calculated based on dE/dX, higher dE/dX assumed to be the end of the track in calorimetry
     // Res range reversed:
     //   1 - residual range calculated from the vertex of the reco track;
     //   0 - residual range calculated traditionally from the end of the reco track;
     int range_reversed = (signal->track_range[shortest_trk_idx][0] >
     			   signal->track_range[shortest_trk_idx][npoints-1]);

     //follow Marshall's convention pg 9,10
     // use res_range-based track length - residual range can be longer that reco track length
     double track_length = signal->track_range[shortest_trk_idx][(npoints-1)*!range_reversed];;
     if (npoints < 4) continue; // have nice templates for longer tracks!
     cuts_stats[kShortTooShort]++;

     // NOTE: track_range is residual range
     //   calorimetry points are stored from the end of the
     //   reconstructed track.
     for( int j=0; j < npoints; ++j){
       double dedx = signal->track_dE_dx[shortest_trk_idx][j];
       double range = signal->track_range[shortest_trk_idx][j];

       if ( (sh_trk_dir == -1 && !range_reversed) ||
	    (sh_trk_dir ==  1 &&  range_reversed) )
	 range = track_length - range;

       h_dEdx     -> Fill( range, dedx );
       h_dEdx_rev -> Fill( track_length - range, dedx );
     }

  }/// All entries end
  cout<<"|"<<endl;

  f->Write();
  f->Close();

  const int width = 25;
  cout<<setw(width)<<"name"<<setw(width-15)<<"stat"
      <<setw(7)<<"Rel%"<<setw(5)<<"Tot%"<<endl;
  for (int i = 0; i < 8; i++) {
    double fraction = 100.*(double)cuts_stats[i]/cuts_stats[0];
    double rel = (i>0)?(100.*(double)cuts_stats[i]/cuts_stats[i-1]):100;
    cout<<setw(width)<<cut_names[i]<<setw(width-15)<<cuts_stats[i]
	<<setw(6)<<setprecision(0)<<fixed<<rel<<"%"
	<<setw(4)<<setprecision(0)<<fixed<<fraction<<"%"
	<<endl;
  }


  return 0;
}
//======================================================================
int main( int argc, char *argv[] ){
  cout<<"**************************** "<<endl;
  cout<<"*    WELCOME TO JAMAICA    * "<<endl;
  cout<<"*      HAVE A NICE DAY     * "<<endl;
  cout<<"**************************** "<<endl;

  if( argc == 1){
    cout<<"Enter a ROOT file to process and output.ROOT file name ..."<<endl;
    return 0;
  }
  if( argc == 3) return Template_PIDloglike_range(argv[1],argv[2]);

  else {
      cout<<"Enter a ROOT file to process and output.ROOT file name ..."<<endl;
      return 0;
  }
}
