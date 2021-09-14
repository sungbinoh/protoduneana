#ifndef EVENTSELECTOR_H
#define EVENTSELECTOR_H

#define NDKAna_cxx
#include "include/NDKAna.h"

class EventSelector {

public:
    enum cut_t {
	kNoCut,
	kTwoTracks,
	kCommonVtx,
	kNoLongProton,
	kNoShortPion,
	kNoPi0,
	kNCaloPoints,
    };

    const int fNcuts = 7;
    const char* cut_names[7] = {
	"No Cut",
	"Two Tracks",
	"Common Vertex",
	"Long is not P",
	"Short is not Pi",
	"No pi0",
	"N Calo pts"
    };

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

    int fRangeReversed[3];

    int fNpointCut;
public:
    EventSelector(const char* fname, const char* templates_sig, const char* templates_bkg, const char* outfname);

private:
    int FindLongestShortest(NDKAna* signal, int& shortest_trk_idx, int& longest_trk_idx);
    int GetDirectionOfTracks(NDKAna* signal,
			     int shortest_trk_idx, int longest_trk_idx,
			     int& sh_trk_dir, int& lo_trk_dir,
			     int& common_vtx);
    void printStats (unsigned* stats);

public:
    // these to be implemented!!!
    virtual int Init();
    virtual int ProcessEvent(NDKAna* evt);
    virtual int Finalise();

    // commont loop over events, with commmon selection
    int Loop(unsigned nentries = 0);

    // setters and getters
    int  GetNpointCut() {return fNpointCut;}
    void SetNpointCut(int cut) { fNpointCut = cut;}

};  // EventSelector class

EventSelector::EventSelector(const char* fname,
			     const char* templates_sig, const char* templates_bkg,
			     const char* outfname)
    : ffname(fname),
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
      fNpointCut(0)
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

    TChain *tree = new TChain("NDKAna/Event");
    tree->Add(fname);
    fEvent = new NDKAna(tree);
    fTree = tree;
}

int EventSelector::Loop(unsigned nentries)
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
	//======= Reco
	// must have 2 tracks exactly - easier handling
	if( fEvent->n_reco_tracks != 2 ) continue;
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

	int totalnpoints = 0;
	// check ordering of calo points: 1 - last in track, first in calo points; 0 - oposite
	for (int iplane = 0; iplane < 3; iplane++) {
	    int npoints = fEvent->n_cal_points_byplane[fShortIdx][iplane];
	    totalnpoints += npoints;
	    fRangeReversed[iplane] =
		( fEvent->track_range_byplane[shortest_trk_idx][iplane][0] >
		  fEvent->track_range_byplane[shortest_trk_idx][iplane][npoints-1] );
	}

	if (!have_common_vtx) continue;
	cuts_stats[kCommonVtx]++;



	// //longest track must be a muon
	// //short track cannot be a pion or lepton
	// if( abs(fEvent->track_mcPDG[longest_trk_idx]) != 13 || fEvent->track_mcPDG[shortest_trk_idx] < 321 ) continue;
	//longest track can't be a proton
	if( abs(fEvent->track_mcPDG[longest_trk_idx]) >= 2212 ) continue;
	cuts_stats[kNoLongProton]++;
	//short track cannot be a pion or lepton
	if( abs(fEvent->track_mcPDG[shortest_trk_idx]) < 321 ) continue;
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

	//==== Cut on n calorimetric points
	if (totalnpoints < fNpointCut) continue;
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

int EventSelector::FindLongestShortest(NDKAna* signal, int& shortest_trk_idx, int& longest_trk_idx) {

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
