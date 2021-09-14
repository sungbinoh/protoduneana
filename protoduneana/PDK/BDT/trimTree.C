#include "EventSelector.h"

class NewEventSelector : public EventSelector {

public:
    NewEventSelector(const char* fname, const char* templates_sig, const char* templates_bkg, const char* outfname)
	: EventSelector(fname, templates_sig, templates_bkg, outfname) {}

    virtual int Init();
    virtual int ProcessEvent(NDKAna* evt);
    virtual int Finalise();

public:
    TTree* fNewTree;

};

int trimTree(const char* fname, const char* outfname)
{

    NewEventSelector sel(fname, "", "", outfname);

    sel.Loop();

    return 0;
}

int NewEventSelector::Init()
{
    foutf = new TFile(foutfname,"recreate");
    fNewTree = fTree->CloneTree(0);
    return 0;
}


int NewEventSelector::ProcessEvent(NDKAna* evt)
{
    fNewTree->Fill();
    return 0;
}

int NewEventSelector::Finalise()
{
    cout<< "Storing " << fNewTree->GetEntries()
	<< " entries into the trimmed tree." << endl;

    foutf -> Close();
    return 0;
}

int EventSelector::Init() {return 0;}
int EventSelector::ProcessEvent(NDKAna* evt) {return 0;}
int EventSelector::Finalise() {return 0;}
