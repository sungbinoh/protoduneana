#ifndef ANALYSER_H
#define ANALYSER_H

#define NDKAna_cxx
#include "NDKAna.h"

class Analyser {

public:
    enum cut_t {
	kNoCut,
	kIsKmu2,
	kAtLeastTwoTracks,
	kTwoTracks,
	kCommonVtx,
	kMuLenght,
	kNoLongProton,
	kNoShortPion,
	kNoPi0,
	kNCaloPoints,
    };

    static const int fNcuts = 10;
    const char* cut_names[fNcuts] = {
	"No Cut",
	"Is Kmu2",
	"At least 2 tracks",
	"Two Tracks",
	"Common Vertex",
	"Mu track length",
	"Long is not P",
	"Short is not Pi",
	"No pi0",
	"N Calo pts",
    };

private:
    const char* ffname;
    const char* foutfname;

    TFile* foutf;

    std::map<TString, TH1*> fHists;

    NDKAna* fEvent;

    TTree* fTree;

    int fShortIdx;
    int fLongIdx;
    int fShortDir;
    int fLongDir;

    int fRangeReversed;

    size_t fEvtIdx;

public:
    unsigned fCutMask;

public:
    Analyser(const char* fname, const char* outfname);
    void AllowBranches(const std::vector<string>& br_list);

private:
    int FindLongestShortest(NDKAna* signal, int& shortest_trk_idx, int& longest_trk_idx);
    int GetDirectionOfTracks(NDKAna* signal,
			     int shortest_trk_idx, int longest_trk_idx,
			     int& sh_trk_dir, int& lo_trk_dir,
			     int& common_vtx);
    void printStats (unsigned* stats);

    bool IsKmu2();

public:
    // these to be implemented!!!
    int Init();
    int ProcessEvent(NDKAna* evt);
    int Finalise();

    // commont loop over events, with commmon selection
    int Loop(unsigned nentries = 0);

    NDKAna* GetEvent() {return fEvent;}
    TTree* GetTree() {return fTree;}

};  // Analyser class

Analyser::Analyser(const char* fname,
		   const char* outfname)
    : ffname(fname),
      foutfname(outfname),
      foutf(0),
      fEvent(0),
      fTree(0),
      fShortIdx(-1),
      fLongIdx(-1),
      fShortDir(0),
      fLongDir(0),
      fCutMask(0)
{

    //========================================================
    // input data file(s)
    //========================================================

    // TChain *tree = new TChain("NDKAna/Event");
    // tree->Add(fname);

    fEvent = new NDKAna(fname);
    fTree = fEvent->fChain;
    fTree->SetBranchStatus("*",0);

    // allow only necessary branches
    vector<string> allowed ({"n_reco_tracks"
		, "track_length"
		, "n_cal_points"
		, "track_vtx"
		, "track_end"
		, "mc_npart"
		, "mc_mother"
		, "mc_pdg"
		, "track_mcPDG"
		, "track_range"
		, "track_dE_dx"});
    AllowBranches(allowed);

}

void Analyser::AllowBranches(const std::vector<string>& br_list)
{
    for (auto it = br_list.begin(); it != br_list.end(); it++) {
	fTree -> SetBranchStatus(it->c_str(),1);
    }
    return;
}


int Analyser::Loop(unsigned nentries)
{
    Init();

    //==================================================
    // Cut stats
    unsigned int cuts_stats[10] = {};


    //========================================================
    float n_entries = fTree->GetEntries();
    if (nentries > 0 && nentries < n_entries)
	n_entries = nentries;
    //n_entries = 10000;
    cout<<"..this many entries "<<n_entries<<endl;
    int fiftieth = n_entries / 50;
    cout<<"|                                                  |\r|";
    for( size_t i = 0; i < n_entries; ++i){
	fEvtIdx = i;
	if ( (i+1) % fiftieth == 0 ) {
	    cout<<"-";
	    cout.flush();
	}

	fEvent->GetEntry(i);

	int longest_trk_idx =-1;
	int shortest_trk_idx =-1;
	FindLongestShortest(fEvent, shortest_trk_idx, longest_trk_idx);
	fShortIdx = shortest_trk_idx;
	fLongIdx = longest_trk_idx;

	cuts_stats[kNoCut]++;

	if ((fCutMask & (1<<kIsKmu2)) && !IsKmu2()) continue;
	cuts_stats[kIsKmu2]++;


	//======= Reco
	// must have at least 2 tracks
	if( (fCutMask & (1<<kAtLeastTwoTracks)) && fEvent->n_reco_tracks < 2 ) continue;
	cuts_stats[kAtLeastTwoTracks]++;

	// must have 2 tracks exactly - easier handling
	if( (fCutMask & (1<<kTwoTracks)) && fEvent->n_reco_tracks != 2 ) continue;
	cuts_stats[kTwoTracks]++;

	int sh_trk_dir =0;
	int lo_trk_dir =0;
	int have_common_vtx = 0;
	GetDirectionOfTracks(fEvent,
			     shortest_trk_idx, longest_trk_idx,
			     sh_trk_dir, lo_trk_dir,
			     have_common_vtx);
	fShortDir = sh_trk_dir;
	fLongDir = lo_trk_dir;

	if ((fCutMask & (1<<kCommonVtx)) && !have_common_vtx) continue;
	cuts_stats[kCommonVtx]++;


	// cut on length of the longest track
	double longest_length = fEvent->track_length[fLongIdx];
	// cut based on Chris Marshall's definition from chris_20180927_pknuLL.pdf
	if ((fCutMask & (1<<kMuLenght)) && (longest_length < 45 || longest_length > 65)) continue;
	cuts_stats[kMuLenght]++;


	// //longest track must be a muon
	// //short track cannot be a pion or lepton
	// if( abs(fEvent->track_mcPDG[longest_trk_idx]) != 13 || fEvent->track_mcPDG[shortest_trk_idx] < 321 ) continue;
	//longest track can't be a proton
	if(  (fCutMask & (1<<kNoLongProton)) &&  abs(fEvent->track_mcPDG[longest_trk_idx]) == 2212 )
	    continue;
	cuts_stats[kNoLongProton]++;
	//short track cannot be a pion or lepton
	if( (fCutMask & (1<<kNoShortPion)) && abs(fEvent->track_mcPDG[shortest_trk_idx]) < 321 )
	    continue;
	cuts_stats[kNoShortPion]++;



	//==== Cut on true pi0
	int has_pi0 = 0;
	for (int j = 0; j < fEvent->mc_npart; j++) {
	    if (fEvent->mc_pdg[j] == 111) {// this is pi0
		has_pi0 = 1;
		break;
	    }
	}
	cuts_stats[kNoPi0]++;

	int npoints = fEvent->n_cal_points[shortest_trk_idx];

	// check ordering of calo points: 1 - last in track, first in calo points; 0 - oposite
	if (npoints > 1)
	    fRangeReversed = (fEvent->track_range[shortest_trk_idx][0] <
			      fEvent->track_range[shortest_trk_idx][npoints-1]);
	else
	    fRangeReversed = -1;

	//==== Cut on n calorimetric points
	if ((fCutMask & (1<<kNCaloPoints)) && npoints < 5) continue;
	cuts_stats[kNCaloPoints]++;

	// End of event selection


	// Process the event
	ProcessEvent(fEvent);
    }
    cout<<"|"<<endl;

    Finalise();

    printStats(cuts_stats);

    return 0;
}

int Analyser::FindLongestShortest(NDKAna* signal, int& shortest_trk_idx, int& longest_trk_idx) {

    //================================================
    // longest and shortest tracks
    //================================================
    float length_tmp = -999.0;
    float short_tmp = 999999;
    int longest_vtx_ID =-999;
    for( size_t ii=0; ii<signal->n_reco_tracks; ++ii){
        if( signal->track_length[ii] > length_tmp ){
	    length_tmp = signal->track_length[ii];
	    longest_trk_idx = ii;
        }
        if( signal->track_length[ii] < short_tmp){
	    short_tmp = signal->track_length[ii];
	    shortest_trk_idx = ii;
        }
    }

    return 0;
}

int Analyser::GetDirectionOfTracks(NDKAna* signal,
					 int shortest_trk_idx, int longest_trk_idx,
					 int& sh_trk_dir, int& lo_trk_dir,
					 int& common_vtx)
{
    //==================================================
    // find direction of the track
    //==================================================

    TVector3 sh_trk_vtx(signal->track_vtx[shortest_trk_idx][0],signal->track_vtx[shortest_trk_idx][1],signal->track_vtx[shortest_trk_idx][2]);
    TVector3 sh_trk_end(signal->track_end[shortest_trk_idx][0],signal->track_end[shortest_trk_idx][1],signal->track_end[shortest_trk_idx][2]);
    TVector3 lo_trk_vtx(signal->track_vtx[longest_trk_idx][0],signal->track_vtx[longest_trk_idx][1],signal->track_vtx[longest_trk_idx][2]);
    TVector3 lo_trk_end(signal->track_end[longest_trk_idx][0],signal->track_end[longest_trk_idx][1],signal->track_end[longest_trk_idx][2]);

    // how far is long track's vertex from the short track? if the
    // end point is closer, than assume wrong rec direction of long
    // track
    double min_start = std::min((sh_trk_vtx-lo_trk_vtx).Mag(), (sh_trk_end-lo_trk_vtx).Mag() );
    double min_end = std::min((sh_trk_vtx-lo_trk_end).Mag(), (sh_trk_end-lo_trk_end).Mag() );

    // if closest ends within 2 cm, then tag as common vertex
    if (std::min(min_start, min_end) > 5.) common_vtx = 0;
    else common_vtx = 1;

    //check muon and short track directions & common vertex?
    // dir = 1 correct
    // dir =-1 reverse dir
    if( min_start < min_end ){
	// proper long track direction reconstructed
	lo_trk_dir =1;
	// check short track's end point
	if( (sh_trk_vtx-lo_trk_vtx).Mag() < (sh_trk_end - lo_trk_vtx).Mag() ) sh_trk_dir = -1;
	else sh_trk_dir = 1;
    }
    else{
	// wrong long track direction reconstructed
	lo_trk_dir =-1;
	if( (sh_trk_vtx-lo_trk_end).Mag() < (sh_trk_end - lo_trk_end).Mag() ) sh_trk_dir = -1;
	else sh_trk_dir = 1;
    }

    return 0;
}

void Analyser::printStats (unsigned* stats)
{
    // const int width = 15;
    // cout<<setw(width)<<"name"<<setw(width)<<"stat"<<endl;
    // for (int i = 0; i < fNcuts; i++) {
    // 	cout<<setw(width)<<cut_names[i]<<setw(width)<<stats[i]<<endl;
    // }
    const int width = 25;
    cout<<setw(width)<<"name"<<setw(width-15)<<"stat"
	<<setw(7)<<"Rel%"<<setw(5)<<"Tot%"<<endl;
    for (int i = 0; i < fNcuts; i++) {
	double fraction = 100.*(double)stats[i]/stats[0];
	double rel = (i>0)?(100.*(double)stats[i]/stats[i-1]):100;

	cout<<setw(width)<<cut_names[i]<<setw(width-15)<<stats[i]
	    <<setw(6)<<setprecision(0)<<fixed<<rel<<"%"
	    <<setw(4)<<setprecision(0)<<fixed<<fraction<<"%"
	    <<endl;
    }
}

bool Analyser::IsKmu2() {
    //========== find if this is Kmu2 ==========
    int Kid = -10;
    int Kidx = -1;
    int isMu = 0;
    int isNu = 0;
    // find track id of the kaon
    for (int i = 0; i < fEvent->mc_npart; i++) {
	if (fEvent->mc_pdg[i] == 321) {
	    Kid = fEvent->mc_id[i];
	    Kidx = i;
	    break;
	}
    }
    if (Kidx < 0) //did not find Kaon (bkg sample)
	return false;

    for (int i = 0; i < fEvent->mc_npart; i++) {
	if (!isMu && fEvent->mc_mother[i] == Kid && fEvent->mc_pdg[i] == -13) {
	    isMu = 1;
	}
	if (!isNu && fEvent->mc_mother[i] == Kid && fEvent->mc_pdg[i] == 14) {
	    isNu = 1;
	}

	if (isNu && isMu) break;
    }

    return isMu && isNu;
    //========================================
}

void NDKAna::Loop()
{}

#endif // ANALYSER_H
