void make_plot_of_low_tension_wires(){

  std::vector< std::string > APA = {
    "UK001",
    "UK002",
    "US001",
    "US002",
    "US003",
    "US004"
  };

  std::vector< std::string > layers = {
    "ULAYER",
    "VLAYER",
    "XLAYER"
  };

  std::vector< std::string > side = {
    "a",
    "b"
  };

  TFile* _file0 = new TFile("../data/tension_measurements.root", "read");

  for (int iapa = 0; iapa < APA.size(); iapa++){

    for (int ilayer = 0; ilayer < layers.size(); ilayer++){

      std::string treeName = "Tension_ProtoDUNE_" + APA.at(iapa) + "_" + layers.at(ilayer);

      TTree* t = (TTree*)_file0->Get(treeName.c_str());

      for (int iside = 0; iside < side.size(); iside++){

        float tension;
        float x_start;
        float x_end;
        float y_start;
        float y_end;

        std::string tensionString = "side_" + side.at(iside) + "_final_tension";

        t->SetBranchAddress(tensionString.c_str(), &tension);
        t->SetBranchAddress("x_start", &x_start);
        t->SetBranchAddress("x_end", &x_end);
        t->SetBranchAddress("y_start", &y_start);
        t->SetBranchAddress("y_end", &y_end);

        TCanvas* c1 = new TCanvas("c1", "c1", 500, 1000);

        std::string titleString =  APA.at(iapa) + " " + layers.at(ilayer) + ", side " + side.at(iside) + ";z (mm);y (mm)";

        TH2D* h_background = new TH2D("h_background", titleString.c_str(), 100, -500, 2800, 100, -500, 6487);
        h_background->GetYaxis()->SetTitleOffset(1.6);
        h_background->Draw();

        TLine* bottom = new TLine(0,0,2300,0);
        TLine* top    = new TLine(0,5988,2300,5988);
        TLine* left   = new TLine(0,0,0,5988);
        TLine* right  = new TLine(2300,0,2300,5988);

        for (int i = 0; i < t->GetEntries(); i++){
          t->GetEntry(i);

          if (tension <= 0) continue;

          if (tension < 3.5){
            TLine *lineLow = new TLine(y_start, x_start, y_end, x_end);
            lineLow->SetLineColor(kRed+1);
            lineLow->Draw("same");
          }
          else if (tension < 4.0){
            TLine *lineMed = new TLine(y_start, x_start, y_end, x_end);
            lineMed->SetLineColor(kOrange+1);
            lineMed->Draw("same");
          }
          else if (tension < 4.5){
            TLine *lineHigh = new TLine(y_start, x_start, y_end, x_end);
            lineHigh->SetLineColor(kGreen+1);
            lineHigh->Draw("same");
          }

        }
        bottom->Draw();
        top->Draw();
        left->Draw();
        right->Draw();

        std::string printString =  "low_tension_wires_"+APA.at(iapa) + "_" + layers.at(ilayer) + "_side_" + side.at(iside)+".png";
        c1->SaveAs(printString.c_str());

        h_background->Delete();

      }
    }
  }
}
