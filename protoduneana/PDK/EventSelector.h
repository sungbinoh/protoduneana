#ifndef EVENTSELECTOR_H
#define EVENTSELECTOR_H

#include <iostream>
#include <iomanip>

#include "TMath.h"
#include "TVector3.h"
#include "TH2F.h"
#include <TRegexp.h>
#include "TROOT.h"
#define NDKAna_cxx
#include "NDKAna.h"

struct pair_info {
    int longer;
    int shorter;
    double long_len;
    double short_len;
    int outwards;
    int vertex_location;

    bool operator()(pair_info& a, pair_info& b) { return a.long_len < b.long_len; }
};


enum {
    // location of a vertex of two matched tracks:
    //   1st long track, 2nd short track
    kVtxVtx,
    kVtxEnd,
    kEndVtx,
    kEndEnd};

using pair_list = std::vector< pair_info >; // list of indices of matched tracks

class EventSelector {
    /**
     * Selector class to do consistent event selection. This one
     * accepts events with at least 2 reconstructed tracks. Search for
     * K-mu or K-pi pairs is done. Based on requirement of common
     * vertex and assumption that muon track is the longest and should
     * point away from the vertex.
     *
     **/

public:
    enum cut_t {
	kNoCut,
	kIsKmu2,
	kIsKpi2,
	kTwoTracks,
	kCommonVtx,
	kNoLongProton,
	kNoShortPion,
	kNoPi0,
	kNCaloPoints,
    };

    static const int fNcuts = 9;
    const char* cut_names[fNcuts] = {
	"No Cut",
	"Is Kmu2",
	"Is Kpi2",
	"Two Tracks",
	"Common Vertex",
	"Long is not P",
	"Short is not Pi",
	"No pi0",
	"N Calo pts"
    };
    unsigned fCutMask;

protected:
    const char* ffname;
    const char* ftempl_bkg;
    const char* ftempl_sig;
    const char* foutfname;

    TFile* foutf;

    std::map<TString, TH1*> fHists;

    NDKAna* fEvent;

    // probability hists
    TH2* fHProbP;
    TH2* fHProbK;
    TH2* fHProbP_rev;
    TH2* fHProbK_rev;

    TTree* fTree;

    int fShortIdx;
    int fLongIdx;
    int fShortDir;
    int fLongDir;

    int fRangeReversed;

    int fNpointCut;
    double fVertexDistCut;

    int fBestPlaneShort;
    int fBestPlaneLong;
public:
    EventSelector(const char* fname, const char* templates_sig, const char* templates_bkg, const char* outfname);

private:
    int GetBestPlane(int idx, int& bestplane, int& npoints);
    int FindLongestShortest(NDKAna* signal, int& shortest_trk_idx, int& longest_trk_idx);
    int GetDirectionOfTracks(NDKAna* signal,
			     int shortest_trk_idx, int longest_trk_idx,
			     int& sh_trk_dir, int& lo_trk_dir,
			     int& common_vtx);
    pair_list GetTrackPairList();
    bool HaveCommonVertex(pair_info& info);
    void PointsOutwards(pair_info& info);

    void printStats (unsigned* stats);
    bool IsKmu2();
    bool IsKpi2();

public:
    // these to be implemented!!!
    virtual int Init();
    virtual int ProcessEvent(NDKAna* evt);
    virtual int Finalise();

    void AllowBranches(const vector<string>& branches);

    // commont loop over events, with commmon selection
    int Loop(unsigned nentries = 0);

    // setters and getters
    int  GetNpointCut() {return fNpointCut;}
    void SetNpointCut(int cut) { fNpointCut = cut;}
    int  GetVertexDistCut() {return fVertexDistCut;}
    void SetVertexDistCut(int cut) { fVertexDistCut = cut;}

};  // EventSelector class

EventSelector::EventSelector(const char* fname, const char* templates_sig, const char* templates_bkg, const char* outfname)  
  : 
      fCutMask(0),
      ffname(fname),
      ftempl_bkg(templates_bkg),
      ftempl_sig(templates_sig),
      foutfname(outfname),
      foutf(0),
      fEvent(0),
      fTree(0),
      fShortIdx(-1),
      fLongIdx(-1),
      fShortDir(0),
      fLongDir(0),
      fNpointCut(0),
      fVertexDistCut(5.)
      //fCutMask(0)
{
    if ( strcmp(templates_sig, "") ) {
	//======= load loglikelihood templates
	TFile *f_bkgd_template = new TFile(templates_bkg,"READ");
	TH2D *h_dEdx_p =  (TH2D*)f_bkgd_template->Get("h_dEdx");
	TH2D *h_dEdx_rev_p =  (TH2D*)f_bkgd_template->Get("h_dEdx_rev");
	h_dEdx_p->SetName("h_dEdx_bkg");
	h_dEdx_rev_p->SetName("h_dEdx_rev_bkg");


	TFile *f_signal_template = new TFile(templates_sig,"READ");
	TH2D *h_dEdx_k =  (TH2D*)f_signal_template->Get("h_dEdx");
	TH2D *h_dEdx_rev_k =  (TH2D*)f_signal_template->Get("h_dEdx_rev");
	h_dEdx_k->SetName("h_dEdx_sig");
	h_dEdx_rev_k->SetName("h_dEdx_rev_sig");


	// Create prob. hists
	int nbins_r    = h_dEdx_p -> GetXaxis()->GetNbins();
	int nbins_dedx = h_dEdx_p -> GetYaxis()->GetNbins();

	fHProbP = (TH2*)h_dEdx_p -> Clone("hProbP");
	fHProbK = (TH2*)h_dEdx_k -> Clone("hProbK");
	fHProbP_rev = (TH2*)h_dEdx_rev_p -> Clone("hProbP_rev");
	fHProbK_rev = (TH2*)h_dEdx_rev_k -> Clone("hProbK_rev");
	TH1* fHProbP_px = (TH1*)h_dEdx_p->ProjectionX("hProbP_px", 1, nbins_dedx);
	TH1* fHProbK_px = (TH1*)h_dEdx_k->ProjectionX("hProbK_px", 1, nbins_dedx);
	TH1* fHProbK_rev_px = (TH1*)h_dEdx_rev_k->ProjectionX("hProbK_rev_px", 1, nbins_dedx);
	TH1* fHProbP_rev_px = (TH1*)h_dEdx_rev_p->ProjectionX("hProbP_rev_px", 1, nbins_dedx);


	const double PROB_LIMIT = 1e-5;
	for (int ir = 1; ir < (nbins_r+1); ir++) {
	    for (int ie = 1; ie < (nbins_dedx+1); ie++) {
		// proton-like
		double prob = h_dEdx_p->GetBinContent(ir,ie)/ fHProbP_px->GetBinContent(ir);
		if (prob < PROB_LIMIT) prob = PROB_LIMIT;
		fHProbP -> SetBinContent(ir, ie, TMath::Log(prob));
		// k-like
		prob = h_dEdx_k->GetBinContent(ir,ie)/ fHProbK_px->GetBinContent(ir);
		if (prob < PROB_LIMIT) prob = PROB_LIMIT;
		fHProbK -> SetBinContent(ir, ie, TMath::Log(prob));
		// reversed
		// proton-like
		prob = h_dEdx_rev_p->GetBinContent(ir,ie)/ fHProbP_rev_px->GetBinContent(ir);
		if (prob < PROB_LIMIT) prob = PROB_LIMIT;
		fHProbP_rev -> SetBinContent(ir, ie, TMath::Log(prob));
		// k-like
		prob = h_dEdx_rev_k->GetBinContent(ir,ie)/ fHProbK_rev_px->GetBinContent(ir);
		if (prob < PROB_LIMIT) prob = PROB_LIMIT;
		fHProbK_rev -> SetBinContent(ir, ie, TMath::Log(prob));
	    }
	}
    }

    //========================================================
    // input data file(s)
    //========================================================

    // TChain *tree = new TChain("NDKAna/Event");
    // tree->Add(fname);
    fEvent = new NDKAna(fname);
    fTree = fEvent->fChain;

    // allow only necessary branches
    vector<string> allowed ({ "subRunNo"
		, "eventNo"
		, "n_reco_tracks"
		, "track_length"
		, "track_vtx"
		, "track_end"
		, "mc_npart"
		, "mc_mother"
		, "mc_pdg"
		, "mc_id"
		, "track_mcPDG"
		, "n_cal_points"
		, "track_range"
		, "track_dE_dx"
		, "n_cal_points_byplane"
		, "track_range_byplane"
		, "track_dE_dx_byplane"});
    fTree->SetBranchStatus("*", 0);
    AllowBranches(allowed);
}

int EventSelector::Loop(unsigned nentries)
{
    Init();
    // print which cuts are turned on
    cout<<"Cuts switched on: ";
    for (int i = 0; i < fNcuts; i++) {
	if ( fCutMask & (1<<i) )
	    cout<<cut_names[i]<<", ";
    }
    cout<<endl;

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
    cout.flush();
    for( int i = 0; i < n_entries; ++i){

	if (fiftieth && (i+1) % fiftieth == 0 ) {
	    cout<<"-";
	    cout.flush();
	}
	//cout<<"Getting entry "<<i<<endl;
	fEvent->GetEntry(i);

	cuts_stats[kNoCut]++;


	// select only Kmu2 events
	if ((fCutMask & (1<<kIsKmu2)) && !IsKmu2()) continue;
	cuts_stats[kIsKmu2]++;

	// select only Kmu2 events
	if ((fCutMask & (1<<kIsKpi2)) && !IsKpi2()) continue;
	cuts_stats[kIsKpi2]++;


	//======= Reco
	// must have at least 2 tracks
	if( (fCutMask & (1<<kTwoTracks)) && fEvent->n_reco_tracks < 2 ) continue;
	cuts_stats[kTwoTracks]++;

	//========================================
	// find matching track pair
	pair_list trk_pairs = GetTrackPairList();


	if ((fCutMask & (1<<kCommonVtx)) && !trk_pairs.size()) continue;
	cuts_stats[kCommonVtx]++;

	static std::map<int, const char*> labels = {{kVtxVtx,"kVtxVtx"},
						    {kVtxEnd,"kVtxEnd"},
						    {kEndVtx,"kEndVtx"},
						    {kEndEnd,"kEndEnd"}, };

	// order track pairs according to the longest track's length
	pair_list trk_pairs_srtd = trk_pairs;
	std::sort(trk_pairs_srtd.begin(), trk_pairs_srtd.end(), trk_pairs.front());

	// pick the pair with the longest long track
	// see if we have multiple pairs with the longest track
	std::vector<int> longest_indices;
	double max = trk_pairs_srtd.back().long_len;
	for (int ii = trk_pairs_srtd.size()-1; ii >= 0; ii--){
	    if (trk_pairs_srtd[ii].long_len < max) break;
	    longest_indices.push_back(ii);
	}
	pair_info selected = trk_pairs_srtd[longest_indices[0]];


	if (longest_indices.size() > 1) {
	    // we have more than 1 candidate! deal with it.
	    // cout<<"More than 1 matched pairs share the longest track. ";
	    // cout<<"SubRun: "<<fEvent->subRunNo
	    // 	<<", Event: "<<fEvent->eventNo
	    // 	<<", Entry: "<<i<<endl;
	    std::vector<pair_info> to_process;
	    for (auto it = longest_indices.begin(); it != longest_indices.end(); it++){
		if (trk_pairs_srtd[*it].outwards)
		    to_process.push_back(trk_pairs_srtd[*it]);
	    }
	    // if (to_process.size() > 1) {
	    // 	cout<<"And "<<to_process.size()<<" pairs have long track pointing outwards. "<<endl;
	    // }
	    if (to_process.size()) {
		// select the one with longest shorter track
		double max_short = to_process[0].short_len;
		selected = to_process[0];
		for (auto p: to_process) {
		    if ( p.short_len > max_short )
			selected = p;
		}
		// cout<<"Selected pair ("<<selected.longer<<","<<selected.shorter<<")"<<endl;
	    }
	}


	fShortIdx = selected.shorter;
	fLongIdx = selected.longer;
	int npoints_short = 0;
	int npoints_long = 0;
	GetBestPlane(fLongIdx, fBestPlaneLong, npoints_long);
	GetBestPlane(fShortIdx, fBestPlaneShort, npoints_short);


	fShortDir = ( selected.vertex_location == kVtxEnd ||
		      selected.vertex_location == kEndEnd ) ? 1: -1;
	fLongDir = ( selected.vertex_location == kVtxEnd ||
		      selected.vertex_location == kVtxVtx ) ? 1: -1;
	//========================================

	// int sh_trk_dir =0;
	// int lo_trk_dir =0;
	// int have_common_vtx = 0;
	// GetDirectionOfTracks(fEvent,
	// 		     shortest_trk_idx, longest_trk_idx,
	// 		     sh_trk_dir, lo_trk_dir,
	// 		     have_common_vtx);



	// check ordering of calo points: 1 - last in track, first in calo points; 0 - oposite
	//cout<<"ShortIDX: "<<fShortIdx
	//    <<", BestPlaneShort: "<<fBestPlaneShort
	//    <<", npoints short: "<<npoints_short<<endl;
	fRangeReversed = (fEvent->track_range_byplane[fShortIdx][fBestPlaneShort][0] >
			  fEvent->track_range_byplane[fShortIdx][fBestPlaneShort][npoints_short-1]);



	// //longest track must be a muon
	// //short track cannot be a pion or lepton
	// if( abs(fEvent->track_mcPDG[fLongIdx]) != 13 || fEvent->track_mcPDG[fShortIdx] < 321 ) continue;
	//longest track can't be a proton
	if( (fCutMask & (1<<kNoLongProton)) && abs(fEvent->track_mcPDG[fLongIdx]) >= 2212 ) continue;
	cuts_stats[kNoLongProton]++;
	//short track cannot be a pion or lepton
	if( (fCutMask & (1<<kNoShortPion)) && abs(fEvent->track_mcPDG[fShortIdx]) < 321 ) continue;
	cuts_stats[kNoShortPion]++;


	//==== Cut on true pi0
	if (fCutMask & (1<<kNoPi0)) {
	    int has_pi0 = 0;
	
             for (int j = 0; j < fEvent->mc_npart; j++) {
		if (fEvent->mc_pdg[j] == 111) {// this is pi0
		    has_pi0 = 1;
		    break;
		}
	    }
	 std::cout<<"has_pi0"<<has_pi0<<endl;	
	}
	cuts_stats[kNoPi0]++;

	//==== Cut on n calorimetric points
	if ((fCutMask & (1<<kNCaloPoints)) && npoints_short < fNpointCut) continue;
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

void EventSelector::AllowBranches(const vector<string>& branches)
{
  for (auto it = branches.begin(); it != branches.end(); it++)
      fTree->SetBranchStatus(it->c_str(), 1);
}

int EventSelector::GetBestPlane(int idx, int& bestplane, int& npoints)
{
    // find the best plane
    bestplane = 0;
    npoints = fEvent->n_cal_points_byplane[idx][0];
    for (int i = 1; i < 3; i++) {
	if (fEvent->n_cal_points_byplane[idx][i] > npoints) {
	    bestplane = i;
	    npoints = fEvent->n_cal_points_byplane[idx][i];
	}
    }

    return bestplane;
}



int EventSelector::FindLongestShortest(NDKAna* signal, int& shortest_trk_idx, int& longest_trk_idx) {

    //================================================
    // longest and shortest tracks
    //================================================
    float length_tmp = -999.0;
    float short_tmp = 999999;
    int longest_vtx_ID =-999;
    std::cout<<"longest_vtx_ID" <<longest_vtx_ID<<endl;
    for( int ii=0; ii<signal->n_reco_tracks; ++ii){
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

int EventSelector::GetDirectionOfTracks(NDKAna* signal,
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

    // if closest ends within 5 cm, then tag as common vertex
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


pair_list EventSelector::GetTrackPairList()
{
    pair_list result;
    int ntracks = fEvent -> n_reco_tracks;

    for (int i = 0; i < (ntracks-1); i++) {
	for (int j = i+1; j < ntracks; j++) {
	    int longer = i;
	    int shorter = j;
	    if (fEvent->track_length[longer] < fEvent->track_length[shorter]) {
		longer = j;
		shorter = i;
	    }
	    pair_info info = {longer, shorter,
			      fEvent->track_length[longer],
			      fEvent->track_length[shorter],
			      0};
	    // check tracks i and j if they have a common vertex
	    if (HaveCommonVertex(info) ) {
		PointsOutwards(info);
		result.push_back(info);
	    }
	}
    }

    return result;
}

bool EventSelector::HaveCommonVertex(pair_info& apair)
// collect info on matched pairs.
{
    static std::map<int, const char*> labels = {{kVtxVtx,"kVtxVtx"},
						{kVtxEnd,"kVtxEnd"},
						{kEndVtx,"kEndVtx"},
						{kEndEnd,"kEndEnd"}, };

    int i = apair.longer;
    int j = apair.shorter;
    TVector3 vtx1(fEvent->track_vtx[i]);
    TVector3 end1(fEvent->track_end[i]);

    TVector3 vtx2(fEvent->track_vtx[j]);
    TVector3 end2(fEvent->track_end[j]);

    double dists[4] = {(vtx1-vtx2).Mag(), (vtx1-end2).Mag(),
		       (end1-vtx2).Mag(), (end1-end2).Mag()};
    int types[4] = { kVtxVtx, kVtxEnd, kEndVtx, kEndEnd };


    double min_dist = dists[0];
    apair.vertex_location = types[0];
    for (int i = 1; i < 4; i++) {
	if (dists[i] < min_dist) {
	    apair.vertex_location = types[i];
	    min_dist = dists[i];
	}
    }

    // // check the minimum distance of 1st track's vtx to 2nd tracks either end
    // double min_dist_vtx = std::min(dists[0], dists[1]);
    // double min_dist_end = std::min(dists[2], dists[3]);
    // // check the closest distance of the 2 track's ends
    // double min_dist = std::min(min_dist_vtx, min_dist_end);

    // if (min_dist_vtx < min_dist_end)
    // 	if (dists[0] < dists[1])
    // 	    apair.vertex_location = kVtxVtx;
    // 	else
    // 	    apair.vertex_location = kVtxEnd;
    // else
    // 	if (dists[2] < dists[3])
    // 	    apair.vertex_location = kEndVtx;
    // 	else
    // 	    apair.vertex_location = kEndEnd;


    if (min_dist < fVertexDistCut) // tracks' end within predefined distance cut
	return 1;

    return 0;
}

void EventSelector::PointsOutwards(pair_info& apair)
{
    int idx_long = apair.longer;

    double sumdedx_vtx = 0;
    double sumdedx_end = 0;
    int npts = 0;
    int bestplane = 0;
    GetBestPlane(idx_long, bestplane, npts);
    double* calo_pts = fEvent->track_dE_dx_byplane[idx_long][bestplane];
    int howmany = npts/2;
    if (howmany > 6) howmany = 6;

    for (int ii = 0; ii < howmany; ii++) {
	sumdedx_vtx += calo_pts[npts-ii-1]; // last calo points are from the track's vertex
	sumdedx_end += calo_pts[ii];
    }

    // if (fEvent->eventNo == 191 ) {
    // 	cout << "pair ("<< idx_long << "," << apair.shorter << "):"
    // 	     << " long best plane: " << fBestPlaneLong << " npoints: "<< npts
    // 	     << " sum dedx at vtx: " << sumdedx_vtx
    // 	     << " sum dedx at end: " << sumdedx_end;
    // }
    apair.outwards = 1;
    int where_vtx = (apair.vertex_location == kEndEnd || apair.vertex_location == kEndVtx);
    if ( (where_vtx && sumdedx_end > sumdedx_vtx) || (!where_vtx && sumdedx_end < sumdedx_vtx))
	// common vertex at the end of the track and sum of dEdx is larger there, too
	// or common vertex in the vertex of the track and the sumdedx is larger there
	// thus, the track looks oriented towards the common vertex
	apair.outwards = 0;

    // if (fEvent->eventNo == 191 ) {
    // 	cout << " vertex at end: " << where_vtx
    // 	     << " outwards: " << apair.outwards << endl;
    // }
}

bool EventSelector::IsKmu2() {
    //========== find if this is Kmu2 ==========
    vector<int> Kid;
    vector<int> Kidx;

    // find track id of the kaon
    for (int i = 0; i < fEvent->mc_npart; i++) {
	if (fEvent->mc_pdg[i] == 321) {
	    Kid.push_back(fEvent->mc_id[i]);
	    Kidx.push_back(i);
	}
    }

    if (!Kidx.size()) //did not find Kaon (bkg sample)
	return false;


    int iskmu2 = 0;
    for (size_t ik = 0; ik < Kidx.size(); ik++) {
	int id = Kid[ik];
	int idx = Kidx[ik];
	std::cout<<"idx"<<idx<<endl;
        int isMu = 0;
	int isNu = 0;
	int ndaughters = 0;
	for (int i = 0; i < fEvent->mc_npart; i++) {
	    if (fEvent->mc_mother[i] != id) continue;
	    ndaughters++;
	    if (!isMu && fEvent->mc_pdg[i] == -13)
		isMu = 1;
	    if (!isNu && fEvent->mc_pdg[i] == 14)
		isNu = 1;
	}

	if (isNu && isMu && ndaughters == 2) {
	    iskmu2 = 1;
	    break;
	}
    }

    return iskmu2;
    //========================================
}

bool EventSelector::IsKpi2() {
    //========== find if this is Kmu2 ==========
    vector<int> Kid;
    vector<int> Kidx;

    // find track id of the kaon
    for (int i = 0; i < fEvent->mc_npart; i++) {
	if (fEvent->mc_pdg[i] == 321) {
	    Kid.push_back(fEvent->mc_id[i]);
	    Kidx.push_back(i);
	}
    }

    if (!Kidx.size()) //did not find Kaon (bkg sample)
	return false;


    int iskpi2 = 0;
    for (size_t ik = 0; ik < Kidx.size(); ik++) {
	int id = Kid[ik];
	int idx = Kidx[ik];
	std::cout<<"idx"<<idx<< endl;
        int isPi0 = 0;
	int isPi_plus = 0;
	int ndaughters = 0;
	for (int i = 0; i < fEvent->mc_npart; i++) {
	    if (fEvent->mc_mother[i] != id) continue;
	    ndaughters++;
	    if (!isPi0 && fEvent->mc_pdg[i] == 111)
		isPi0 = 1;
	    if (!isPi_plus && fEvent->mc_pdg[i] == 211)
		isPi_plus = 1;
	}

	if (isPi_plus && isPi0 && ndaughters == 2) {
	    iskpi2 = 1;
	    break;
	}
    }

    return iskpi2;
    //========================================
}

void EventSelector::printStats (unsigned* stats)
{
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



void NDKAna::Loop()
{}

#endif // EVENTSELECTOR_H
