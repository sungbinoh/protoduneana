
TGraphErrors* make_graph(TH2D* hist, std::string type){

  std::vector<double> xval;
  std::vector<double> xvalerr;
  std::vector<double> yval;
  std::vector<double> yvalerr;

  for (int i = 1; i < hist->GetNbinsX(); i++){

    if (hist->Integral(i,i+1) == 0) continue;

    if ( type == "rms" ){
      double rms = 0;      
      double nEntries = 0;
      for (int j = 1; j < hist->GetNbinsY(); j++){
        rms += hist->GetBinContent(i,j)*std::pow(hist->GetYaxis()->GetBinCenter(j),2);
        nEntries += hist->GetBinContent(i,j);
      }

      if (nEntries == 0) continue;

      rms /= nEntries;
      rms = std::sqrt(rms);

      yval.push_back(rms);
      yvalerr.push_back(0);
    }

    else if ( type == "mean" ){
      double mean = 0;
      double nEntries = 0;
      for (int j = 1; j < hist->GetNbinsY(); j++){
        mean += (hist->GetBinContent(i,j)*hist->GetYaxis()->GetBinCenter(j));
        nEntries += hist->GetBinContent(i,j);
      }

      if (nEntries == 0) continue;

      mean /= nEntries;

      double meanerr = 0;
      for (int j = 1; j < hist->GetNbinsY(); j++){
        meanerr += (hist->GetBinContent(i,j)*std::pow(hist->GetYaxis()->GetBinCenter(j) - mean, 2));
      }
      meanerr /= nEntries;
      meanerr = std::sqrt(meanerr);
      meanerr = meanerr/sqrt(nEntries);

      yval.push_back(mean);
      yvalerr.push_back(meanerr);

    }

    else if ( type == "median" ){
      double median = 0;
      std::vector<double> entries;
      for (int j = 1; j < hist->GetNbinsY(); j++){
        for (int k = 0; k < hist->GetBinContent(i,j); k++){
          entries.push_back(hist->GetYaxis()->GetBinCenter(j));
        }
      }

      if (entries.size() == 0) continue;

      median = TMath::Median(entries.size(), &entries[0]);

      double medianerr = 0;

      yval.push_back(median);
      yvalerr.push_back(medianerr);

    }

    else throw std::logic_error("invalid type sent to make_graph");

    xval.push_back(hist->GetXaxis()->GetBinCenter(i));
    xvalerr.push_back(hist->GetXaxis()->GetBinWidth(i)/2.);

  }

  //for (int xx = 0; xx < xval.size(); xx++){
  //  std::cout << "-- bin " << xx << std::endl;
  //  std::cout << xval.at(xx) << " +/- " << xvalerr.at(xx) << std::endl;
  //  std::cout << yval.at(xx) << " +/- " << yvalerr.at(xx) << std::endl;
  //}

  TGraphErrors* gr = new TGraphErrors(xval.size(), &xval[0], &yval[0], &xvalerr[0], &yvalerr[0]);

  return gr;

}

void make_tension_var_plots_per_apa(){

  std::vector< std::string > vars = {
    "trackHitCaloChargeDep",
    "trackHitCaloEnergyDep",
    "trackHitRMS",
    "trackHitPeakTime",
    "trackHitPeakAmplitude",
    "trackHitIntegral"
  };

  std::vector< std::string > labels = {
    ";Wire Tension (N); Hit dQ/dx (ADC*/cm)",
    ";Wire Tension (N); Hit dE/dx (MeV/cm)",
    ";Wire Tension (N); Hit RMS",
    ";Wire Tension (N); Hit Peak Time (#mu s)",
    ";Wire Tension (N); Hit Peak Amplitude (ADC*)",
    ";Wire Tension (N); Hit Integral"
  };

  std::vector< std::vector<int> > bins = {
    {10, 3, 8, 5000, 0, 10000},
    {10, 3, 8, 5000, 0, 500},
    {10, 3, 8, 5000, 0, 100},
    {10, 3, 8, 5000, 0, 10000},
    {10, 3, 8, 5000, 0, 500},
    {10, 3, 8, 5000, 0, 5000}
  };

  std::vector< std::string > cutsView = {
    "trackHitView>-1", // all planes
    "trackHitView==0", // u plane
    "trackHitView==1", // v plane
    "trackHitView==2"  // x plane
  };

  std::vector< std::string > cutsTPC = {
    "trackHitAPABuildNumber>-1", // all tpcs 
    "trackHitAPABuildNumber==1", // tpc 1
    "trackHitAPABuildNumber==2", // tpc 2
    "trackHitAPABuildNumber==3", // tpc 3
    "trackHitAPABuildNumber==4", // tpc 4
    "trackHitAPABuildNumber==5", // tpc 5
    "trackHitAPABuildNumber==6"  // tpc 6
  };

  TFile* fOut = new TFile("2dplotsandgraphs.root", "recreate");

  TTree* tree = (TTree*)_file0->Get("tensionanalysis/analysis_tree");

  for (int ivar = 0; ivar < vars.size(); ivar++){

    for (int icv = 0; icv < cutsView.size(); icv++){

      for (int itpc = 0 ; itpc < cutsTPC.size(); itpc++){

        std::string plotName = std::string(vars.at(ivar)+"_"+cutsView.at(icv)+"_"+cutsTPC.at(itpc));

        TH2D* tensionHist = new TH2D(
            plotName.c_str(), 
            labels.at(ivar).c_str(), 
            bins.at(ivar).at(0), bins.at(ivar).at(1), bins.at(ivar).at(2), 
            bins.at(ivar).at(3), bins.at(ivar).at(4), bins.at(ivar).at(5)); 

        tree->Draw((vars.at(ivar)+std::string(":trackHitWireTension >> ")+plotName).c_str(), (cutsView.at(icv)+std::string(" && ")+cutsTPC.at(itpc)).c_str());

        TCanvas *c1 = new TCanvas();
        tensionHist->SetContour(1000);
        tensionHist->Draw("colz");
        c1->SaveAs((plotName+std::string(".png")).c_str());

        TGraphErrors* rmsGraph = make_graph(tensionHist, "rms");       
        rmsGraph->SetName((plotName+std::string("_rms")).c_str());
        rmsGraph->Write();

        TGraphErrors* medianGraph = make_graph(tensionHist, "median");       
        medianGraph->SetName((plotName+std::string("_median")).c_str());
        medianGraph->Write();

        TGraphErrors* meanGraph = make_graph(tensionHist, "mean");       
        meanGraph->SetName((plotName+std::string("_mean")).c_str());
        meanGraph->Write();

        fOut->Write();

        tensionHist->Delete();

      }

    }
  }
}

