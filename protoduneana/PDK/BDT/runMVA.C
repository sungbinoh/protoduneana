#include "BDTTree.h"

int runMVA(const char* sigfname, const char* bkgfname,
	   const int nTrainSig = 10000, const int nTrainBkg = 10000,
	   const char* outpref = "")
{

    /// Declare Factory
    TMVA::Tools::Instance();

    auto sigf = TFile::Open(sigfname);
    auto bkgf = TFile::Open(bkgfname);
    auto outputFile = TFile::Open( Form("%sTMVAOutputCV_train_S%d_B%d.root",
					outpref, nTrainSig, nTrainBkg), "RECREATE" );

    TMVA::Factory factory("TMVAClassification", outputFile,
			  "!V:ROC:!Correlations:!Silent:Color:"
			  "!DrawProgressBar:AnalysisType=Classification" );


    /// Prepare Data Loader
    TMVA::DataLoader loader("dataset");

#undef DOUBLE
#undef INT
#undef OBSERVER_INT
#define DOUBLE(var) loader.AddVariable(#var)
#define INT(var)    loader.AddVariable(#var)
#define SPECTATOR_INT(var) loader.AddSpectator(#var)

    VAR_LIST;


    /// Setup datasets
    TTree *tsig, *tbkg;
    sigf->GetObject("bdtTree", tsig);
    bkgf->GetObject("bdtTree", tbkg);

    TCut cut_sig, cut_bkg;

    loader.AddSignalTree    (tsig, 1.0);   //signal weight  = 1
    loader.AddBackgroundTree(tbkg, 1.0);   //background weight = 1
    loader.PrepareTrainingAndTestTree(cut_sig, cut_bkg,
				      Form("nTrain_Signal=%d:nTrain_Background=%d"
					   ":SplitMode=Random:NormMode=NumEvents:!V",
					   nTrainSig, nTrainBkg));

    /// Book MVA method -- BDT
    //Boosted Decision Trees
    // factory.BookMethod(&loader,TMVA::Types::kBDT, "BDT_400",
    // 		       "!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=2"
    // 		       ":BoostType=AdaBoost:AdaBoostBeta=0.5"
    // 		       ":UseBaggedBoost:BaggedSampleFraction=0.5"
    // 		       ":SeparationType=GiniIndex:nCuts=20" );
    factory.BookMethod(&loader,TMVA::Types::kBDT, "BDT_1000",
		       "!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=2"
		       ":BoostType=AdaBoost:AdaBoostBeta=0.5"
		       ":UseBaggedBoost:BaggedSampleFraction=0.5"
		       ":SeparationType=GiniIndex:nCuts=20" );

    /// Train the BDT
    factory.TrainAllMethods();

    /// Test and Evaluate
    factory.TestAllMethods();
    factory.EvaluateAllMethods();


    /// Plot ROC curve
    auto c1 = factory.GetROCCurve(&loader);
    c1->Draw();

    return 0;
}
