R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

void ReadSim() {

    //Histogram of the angular acceptance
    TH1D* TrueAngAcceptance = new TH1D("TrueAngAcceptance","Angular acceptance against cosine of incident angle",100,0,1);

    for (int n=0; n<=100; n++) {
        //create a string of the file name
        std::ostringstream oss;
        oss << "wcsim";
        oss << std::fixed << std::setprecision(2) << n*0.01;
        oss << ".root";

        //access to the file
        TFile* file = new TFile(oss.str().c_str(), "read");

        //access to the tree and branch
        TTree* tree = (TTree*)file->Get("wcsimT");
        WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
        tree->SetBranchAddress("wcsimrootevent", &wcsimrootsuperevent);
        WCSimRootTrigger* wcsimrootevent = 0;

        int nEntries = tree->GetEntries();
        std::cout << "Reading file: " << oss.str() << "\nNumber of entries: " << nEntries << "\n";

        for (int i=0; i<nEntries; i++) {
            delete wcsimrootsuperevent;
            wcsimrootsuperevent = 0;

            tree->GetEntry(i);
            wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
            int THPMTNum = wcsimrootevent->GetNcherenkovhits();

            for (int j=0; j<THPMTNum; j++) {
                WCSimRootCherenkovHit* WRCH = (WCSimRootCherenkovHit*) wcsimrootevent->GetCherenkovHits()->At(j);
                if (WRCH->GetTubeID()!=1103) {continue;}    //get the target PMT
                int FirstHitPos = WRCH->GetTotalPe(0), LastHitPos = WRCH->GetTotalPe(0)+WRCH->GetTotalPe(1)-1;
                //loop for all true hits of the targeted PMT
                for (int k=FirstHitPos; k<=LastHitPos; k++) {
                    WCSimRootCherenkovHitTime* WRCHT = (WCSimRootCherenkovHitTime*) wcsimrootevent->GetCherenkovHitTimes()->At(k);
                    if (WRCHT->GetParentSavedTrackID()!=0) {continue;}   //include only main photon track
                    TrueAngAcceptance->Fill(TrueAngAcceptance->GetBinCenter(n+1));  //Filling the histogram with the number of true hits being the angular acceptance
                }
            }
        }

        file->Close();
    }

    //create another histogram with the averaged angular acceptance with a bin width of cosine of incident angle = 0.1
    TH1D* AvgAngAcceptance = new TH1D("AvgAngAcceptance","Angular acceptance against cosine of incident angle",10,0,1);
    for (int i=1; i<=100; i++) {
        AvgAngAcceptance->Fill(TrueAngAcceptance->GetBinCenter(i),TrueAngAcceptance->GetBinContent(i));
    }

    //Normalization of the angular acceptance histogram
    double NormalizationFactor = TrueAngAcceptance->GetBinContent(100);
    for (int i=1; i<=100; i++) {
        double CurrentValue = TrueAngAcceptance->GetBinContent(i);
        TrueAngAcceptance->SetBinContent(i, CurrentValue/NormalizationFactor);
    }
    
    //Normalization of the averaged angular acceptance histogram
    NormalizationFactor = AvgAngAcceptance->GetBinContent(10);
    for (int i=1; i<=10; i++) {
        double CurrentValue = AvgAngAcceptance->GetBinContent(i);
        AvgAngAcceptance->SetBinContent(i, CurrentValue/NormalizationFactor);
    }

    //Save the averaged angular acceptance histogram into a root file
    TFile* Hist = new TFile("SimAngAcceptance.root", "RECREATE");
    AvgAngAcceptance->Write();
    Hist->Close();
    
    //plot the angular acceptance histogram
    TCanvas* c1 = new TCanvas();
    TrueAngAcceptance->SetStats(0);
    TrueAngAcceptance->Draw();
    c1->SaveAs(Form("SimTrueAngAcceptance.pdf"));
    delete AvgAngAcceptance;
    delete TrueAngAcceptance;
    delete c1;
}