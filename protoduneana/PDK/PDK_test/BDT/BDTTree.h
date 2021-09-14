#ifndef BDTTREE_H
#define BDTTREE_H

#include "TTree.h"

namespace BDTTree {

    // define all variables to be used here
#define VAR_LIST							\
    INT(n_tracks);							\
    DOUBLE(len_long);							\
    DOUBLE(len_short);							\
    /* DOUBLE(template_dllr); */					\
    INT(short_ncalpts);							\
    DOUBLE(template_sig_forward);					\
    DOUBLE(template_sig_backward);					\
    DOUBLE(template_bkg_forward);					\
    DOUBLE(template_bkg_backward);					\
    DOUBLE(long_pida);							\
    DOUBLE(short_pida);							\
									\
    INT(n_showers);							\
    /*DOUBLE(sh_dist_vtx); // distance of the nearest shower to the p decay vertex */ \
    /*DOUBLE(sh_dist_mu_vtx);  // distance of the nearest shower to the K decay vertex */ \
    /*DOUBLE(sh_dist_mu_end);  // distance of the nearest shower to the muon decay vertex */ \
									\
    DOUBLE(trk_en); /* // distance of the nearest shower to the muon decay vertex */ \
    DOUBLE(sh_en); /* // distance of the nearest shower to the muon decay vertex */ \
    SPECTATOR_INT(subRun);						\
    SPECTATOR_INT(event)							\


#define DOUBLE(var) double var
#define INT(var)    int var
#define SPECTATOR_INT(var) INT(var)

    struct Evt_t {
	VAR_LIST;
    };


#undef DOUBLE
#undef INT
#define DOUBLE(var) tree->Branch(#var, &evt.var, #var"/D")
#define INT(var)    tree->Branch(#var, &evt.var, #var"/I")

    TTree* buildTree(TTree* tree, Evt_t& evt) {
	VAR_LIST;
	return tree;
    }

#undef DOUBLE
#undef INT
#define DOUBLE(var) tree->SetBranchAddress(#var, &evt.var)
#define INT(var)    tree->SetBranchAddress(#var, &evt.var)

    TTree* registerTree(TTree* tree, Evt_t& evt) {
	VAR_LIST;
	return tree;
    }


} // BTDTree namespace
#endif
