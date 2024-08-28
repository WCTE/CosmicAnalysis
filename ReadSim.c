R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

void ReadSim() {

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

            //std::cout<<"THPMTNum = "<<THPMTNum<<std::endl;
            for (int j=0; j<THPMTNum; j++) {
                WCSimRootCherenkovHit* WRCH = (WCSimRootCherenkovHit*) wcsimrootevent->GetCherenkovHits()->At(j);
                if (WRCH->GetTubeID()!=1103) {continue;}
                int FirstHitPos = WRCH->GetTotalPe(0), LastHitPos = WRCH->GetTotalPe(0)+WRCH->GetTotalPe(1)-1;
                for (int k=FirstHitPos; k<=LastHitPos; k++) {
                    WCSimRootCherenkovHitTime* WRCHT = (WCSimRootCherenkovHitTime*) wcsimrootevent->GetCherenkovHitTimes()->At(k);
                    if (WRCHT->GetParentSavedTrackID()!=0) {continue;}   //include only main photon track
                    TrueAngAcceptance->Fill(TrueAngAcceptance->GetBinCenter(n+1));
                }
            }
        }

        file->Close();
    }

    TH1D* AvgAngAcceptance = new TH1D("AvgAngAcceptance","Angular acceptance against cosine of incident angle",10,0,1);
    for (int i=1; i<=100; i++) {
        AvgAngAcceptance->Fill(TrueAngAcceptance->GetBinCenter(i),TrueAngAcceptance->GetBinContent(i));
    }

    double NormalizationFactor = TrueAngAcceptance->GetBinContent(100);
    for (int i=1; i<=100; i++) {
        double CurrentValue = TrueAngAcceptance->GetBinContent(i);
        TrueAngAcceptance->SetBinContent(i, CurrentValue/NormalizationFactor);
    }
    
    NormalizationFactor = AvgAngAcceptance->GetBinContent(10);
    for (int i=1; i<=10; i++) {
        double CurrentValue = AvgAngAcceptance->GetBinContent(i);
        AvgAngAcceptance->SetBinContent(i, CurrentValue/NormalizationFactor);
    }

    TFile* Hist = new TFile("SimAngAcceptance.root", "RECREATE");
    AvgAngAcceptance->Write();
    Hist->Close();
    
    TCanvas* c1 = new TCanvas();
    TrueAngAcceptance->SetStats(0);
    TrueAngAcceptance->Draw();
    c1->SaveAs(Form("SimTrueAngAcceptance.pdf"));
    AvgAngAcceptance->SetStats(0);
    AvgAngAcceptance->Draw("hist");
    c1->SaveAs(Form("AvgAngAcceptance.pdf"));
    delete AvgAngAcceptance;
    delete TrueAngAcceptance;
    delete c1;
}