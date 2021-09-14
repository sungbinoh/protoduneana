//#include "Analyser.h"
//#include "EventDisplay.h"
#include  "NDKAna.h"
//#include  "Processor.h"
#include "EventSelector.h"
#include "BDTTree.h"

class MySelector : public EventSelector
{
public:
    MySelector(const char* fname, const char* templates_sig, const char* templates_bkg, const char* outfname):
	EventSelector(fname, templates_sig, templates_bkg, outfname){}

    virtual int Init();
    virtual int ProcessEvent(NDKAna* evt);
    virtual int Finalise();

    double NearestShower(const TVector3& pos);

    TTree* fOutTree;
    BDTTree::Evt_t fEvt;
};



int createBDTinput(const char* fname, const char* templates_sig,
		   const char* templates_bkg, const char* outfname,
		   int npointcut = 3)
{
    MySelector sel(fname, templates_sig, templates_bkg, outfname);
    sel.SetNpointCut(npointcut);

    sel.Loop();

    return 0;
}


int MySelector::Init()
{
    cout<<"Setting cut mask"<<endl;
    //fCutMask |= (1<<kIsKmu2);
    fCutMask |= (1<<kTwoTracks);
    fCutMask |= (1<<kCommonVtx);
    //fCutMask |= (1<<kNCaloPoints);
    // fCutMask |= (1<<kNoShortPion);

    // allow to read the shower specific branches
    vector<string> allowed ({"n_showers",
		"sh_energy",
		"sh_start_X",
		"sh_start_Y",
		"sh_start_Z",
		"track_PIDA",
		"track_pitch",
		});
    AllowBranches(allowed);

    //==================================================
    // Output file and Hitograms
    //==================================================
    // gStyle->SetPalette(kBlueRedYellow);
    // gStyle->SetNumberContours(255);
    // gStyle->SetHistFillStyle(0);
    // gStyle->SetOptStat(0);

    if ( strcmp(foutfname,"") )
	foutf = new TFile(foutfname,"recreate");

    // create output tree
    fOutTree = new TTree("bdtTree", "BDT tree");
    BDTTree::buildTree(fOutTree, fEvt);

    return 0;
}

int MySelector::ProcessEvent(NDKAna* evt)
{
    //==================================================
    // Calculate log-likelihood from templates
    //==================================================
    double like_p = 0;
    double like_k = 0;
    //double loglike = 0;
    double like_p_rev = 0;
    double like_k_rev = 0;
    //double loglike_rev =0;
    //double PID_loglike = 0;

    double track_length = evt->track_length[fShortIdx];
    std::cout<<"track_lenght"<<track_length<<endl;
    int npoints = evt->n_cal_points_byplane[fShortIdx][fBestPlaneShort];

    // if (npoints < 3 ) return 0;

    static const double max_range = fHProbP -> GetXaxis() -> GetXmax();
    std::cout<<"max_range"<<max_range<<endl;
    double range_track_length = evt->track_range_byplane[fShortIdx][fBestPlaneShort][(npoints-1)*!fRangeReversed];
    int max_r_bin = fHProbP ->GetXaxis()->GetNbins();
    int max_e_bin = fHProbP ->GetYaxis()->GetNbins();


    // compare reconstructed dE/dX vs wire number distance from muon vertex with the templates
    for( int k = 0; k < npoints; k++ ){
	// loop over calorimetric points of the best plane

	// indexing from the muon vertex (0 at the vertex)
	// calo points stored from end of a track - res range = 0 in the first calo point
	int idx = k;


	 double res_range = evt->track_range_byplane[fShortIdx][fBestPlaneShort][idx];
	 if ( (fShortDir == -1 && !fRangeReversed) ||
	      (fShortDir ==  1 && fRangeReversed) )
	     // reco direction swapped (heading out of the vertex),
	     // residual range will be the other end
	     res_range = range_track_length - res_range;//track_length - res_range;
	 if (res_range < 0) res_range = 0;
	 double dEdx =evt->track_dE_dx_byplane[fShortIdx][fBestPlaneShort][idx];

	 int bin_p = fHProbP ->GetYaxis()->FindBin(dEdx);
	 int bin_k = fHProbK ->GetYaxis()->FindBin(dEdx);// required only if tepmlates expected to be binned differently
	 int bin_r = fHProbP ->GetXaxis()->FindBin(res_range);// expect the two templates to have the same range

	 if ( bin_r > 0 && bin_r <= max_r_bin) {
	     double tmp_like_p = fHProbP->GetBinContent(bin_r, bin_p);
	     double tmp_like_k = fHProbK->GetBinContent(bin_r, bin_k);

	     // deal with calo pts where dedx is out of our range - assign
	     // them default minimal likelihood
	     if (bin_p == 0 || bin_p > max_e_bin ) {
		 tmp_like_p = -5;
		 tmp_like_k = -5;
	     }
	     like_p += tmp_like_p;
	     like_k += tmp_like_k;
	 }

	 //using reverse template
	 // int rev_idx = npoints - idx - 1;
	 // double dEdx_rev = evt->track_dE_dx[fShortIdx][rev_idx];
	 double res_range_rev = range_track_length - res_range; //= evt->track_range[fShortIdx][rev_idx];
	 // if (fShortDir == 1) // translate res range to res range from the other end
	 //     res_range_rev = track_length - res_range_rev;
	 // if (res_range_rev < 0) res_range_rev = 0;

	 int bin_p_rev = fHProbP_rev->GetYaxis()->FindBin(dEdx); // required only if tepmlates expected to be binned differently
	 int bin_k_rev = fHProbK_rev->GetYaxis()->FindBin(dEdx); // dtto
	 int bin_r_rev = fHProbP_rev->GetXaxis()->FindBin(res_range_rev);// expect the two templates to have the same range

	 if ( bin_r_rev > 0 && bin_r_rev <= max_r_bin) {
	     double tmp_like_p_rev = fHProbP_rev->GetBinContent(bin_r_rev, bin_p_rev);
	     double tmp_like_k_rev = fHProbK_rev->GetBinContent(bin_r_rev, bin_k_rev);

	     if (bin_p_rev == 0 || bin_p_rev > max_e_bin ) {
		 tmp_like_p_rev = -5;
		 tmp_like_k_rev = -5;
	     }
	     like_p_rev += tmp_like_p_rev;
	     like_k_rev += tmp_like_k_rev;
	 }

    }

    //loglike = like_p-like_k;
   // loglike_rev = like_p_rev-like_k_rev;
  //  PID_loglike = loglike + loglike_rev;

	//std::cout<<"PID_loglike"<<PID_loglike<<endl;

    fEvt.short_ncalpts = npoints;
    //fEvt.template_dllr = PID_loglike;
    fEvt.template_sig_forward = like_k;
    fEvt.template_sig_backward = like_k_rev;
    fEvt.template_bkg_forward = like_p;
    fEvt.template_bkg_backward = like_p_rev;

    fEvt.n_tracks = evt->n_reco_tracks;
    fEvt.len_long = evt->track_length[fLongIdx];
    fEvt.len_short = evt->track_length[fShortIdx];

    // save tracks' PIDA
    fEvt.long_pida = fEvent->track_PIDA[fLongIdx][fBestPlaneLong];
    fEvt.short_pida = fEvent->track_PIDA[fShortIdx][fBestPlaneShort];
    // prevent extremely low PIDA (set to -999 if undefined)
    if (fEvt.long_pida < 0 || isnan(fEvt.long_pida))
	fEvt.long_pida = -1;
    if (fEvt.short_pida < 0 || isnan(fEvt.short_pida))
	fEvt.short_pida = -1;

    fEvt.n_showers = evt->n_showers;

    /*
    // find distances of nearest shower(s) to: 1) proton decay vertex; 2) kaon decay vtx; 3) muon decay vertex
    TVector3 short_vtx(fEvent->track_vtx[fShortIdx]);
    TVector3 short_end(fEvent->track_end[fShortIdx]);
    TVector3 long_vtx(fEvent->track_vtx[fLongIdx]);
    TVector3 long_end(fEvent->track_end[fLongIdx]);

    double dist_vtx = NearestShower(short_vtx);
    double dist_end = NearestShower(short_end);

    fEvt.sh_dist_mu_end = NearestShower((fLongDir==1)?long_end:long_vtx);
    fEvt.sh_dist_vtx = (fShortDir==1)?dist_vtx:dist_end;
    fEvt.sh_dist_mu_vtx = (fShortDir==1)?dist_end:dist_vtx;
    */

    // Record total energy reconstructed in tracks and separately in showers

    // All track energy 1st
    // will be determined from dE/dx of the best plane. Here the best plane
    // is inherited from the original code, and corresponds to the best chi2
    // from the calorimetry module.
    double total_trk_energy = 0.;
    // loop over all Tracks
    for (int itrk = 0; itrk < fEvent->n_reco_tracks; itrk++) {
	// loop over all calo points in the best plane
	for (int icalo = 0; icalo < fEvent->n_cal_points[itrk]; icalo++) {
	    double energy = fEvent->track_dE_dx[itrk][icalo];
	    if (energy > 100.) energy = 100.;
	    energy *= fEvent->track_pitch[itrk][icalo];
	    total_trk_energy += energy;
	}
    }
    fEvt.trk_en = total_trk_energy;

    // All shower energy
    // There is energy stored for each shower, calculated in each plane.
    // Will take the calculation from the best plane.
    double total_shr_energy = 0.;
    for (int ish = 0; ish < fEvent->n_showers; ish++) {
	total_shr_energy += fEvent->sh_energy[ish][fEvent->sh_bestplane[ish]];
    }
    fEvt.sh_en = total_shr_energy;

    // add event identification
    fEvt.subRun = fEvent->subRunNo;
    fEvt.event = fEvent->eventNo;


    fOutTree->Fill();

    return 0;
}


int MySelector::Finalise()
{

    foutf->Write();
    foutf->Close();

    return 0;
}

double MySelector::NearestShower(const TVector3& pos) {
    int size = fEvent->n_showers;

    double min_dist = 9999.;
    for (int i = 0; i<size; i++) {
	TVector3 shower_start (fEvent->sh_start_X[i], fEvent->sh_start_Y[i], fEvent->sh_start_Z[i]);
	double dist = (shower_start-pos).Mag();

	if (dist < min_dist)
	    min_dist = dist;
    }

    if (min_dist == 9999.)
	return -1;

    return min_dist;
}


int EventSelector::Init() {return 0;}
int EventSelector::ProcessEvent(NDKAna* evt) {return 0;}
int EventSelector::Finalise() {return 0;}
