
//int a = 0;
//int subdet =  37;
//int xhn = -99;

{
TFile *fin = TFOile::Open("/dune/app/users/tstokes/v09_22_12/srcs/protoduneana/protoduneana/PDK/ndk_test.root");
fin->ls();
if (fin == 0) {
	printf("Error: cannot open the file!\n");
} else {
	TTree *t = 0;
	fin->GetObject("lemma",t);

	TCanvas *c = new TCanvas("c","canvas", 1280, 1024);
	
	TH1F *TrackLengths = new TH1F("TrackLengths", "", 100, 0., 0., 100, 0., 0.);
	//TString varexp = TString::Format("
	TH1F *track_length_Muons = new TH1F("track_length_Muons", "", 100, 0., 0., 100, 0., 0.);
	TH1F *track_length_Kaons = new TH1F("track_length_Kaons", "", 100, 0., 0., 100, 0., 0.);
	TH1F *track_length_Pions = new TH1F("track_length_Pions", "", 100, 0., 0., 100, 0., 0.);
	//TH1F *track_length_Proton = new TH1F("track_length_Muons", "", 100, 0., 0., 100, 0., 0.);
	t->Draw(track_length_Muons);
	t->Draw(track_length_Kaons);
	t->Draw(track_length_Pions);

c->Print("/dune/app/users/tstokes/v09_22_02/srcs/protoduneana/protoduneana/PDK/EXAMPLE.png");

}
}

