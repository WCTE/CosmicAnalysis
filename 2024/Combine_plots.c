R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

void Combine_plots() {
    TFile* file1 = TFile::Open("./AngAcceptance.root");     //File containing the angular acceptance histogram of the cosmic muon simulation
    TFile* file2 = TFile::Open("./SimAngAcceptance.root");      //File containing the angular acceptance histogram by direct simulation

    TH1D* WCSimTH = (TH1D*)file1->Get("WCSimTH");       //Angular acceptance histogram of WCSim with true hit
    TH1D* fitQunDH = (TH1D*)file1->Get("fitQunDH");     //Angular acceptance histogram of fitQun with digit hit
    TH1D* fitQunTH = (TH1D*)file1->Get("fitQunTH");     //Angular acceptance histogram of fitQun with true hit
    TH1D* AvgAngAcceptance = (TH1D*)file2->Get("AvgAngAcceptance");     //Angular acceptance histogram by direct simulation

    //comment/uncomment the following code to select the histogram to be plotted together
    TCanvas* c1 = new TCanvas();
    WCSimTH->SetStats(0);
    WCSimTH->SetLineColor(kRed);
    fitQunDH->SetStats(0);
    fitQunDH->SetLineColor(kBlue);
    fitQunTH->SetStats(0);
    fitQunTH->SetLineColor(kGreen);
    AvgAngAcceptance->SetStats(0);
    AvgAngAcceptance->SetLineColor(kBlack);
    WCSimTH->SetTitle("Angular acceptance against cosine of incident angle");
    WCSimTH->SetXTitle("cosine");
    WCSimTH->SetYTitle("angular acceptance");
    WCSimTH->Draw("hist");
    fitQunDH->Draw("hist same");
    fitQunTH->Draw("hist same");
    AvgAngAcceptance->Draw("hist same");
    
    TLegend* legend = new TLegend(0.1,0.75,0.48,0.9);
    legend->AddEntry(WCSimTH, "WCSim with true hit", "l");
    legend->AddEntry(fitQunDH, "fitQun with digit hit", "l");
    legend->AddEntry(fitQunTH, "fitQun with true hit", "l");
    legend->AddEntry(AvgAngAcceptance, "Direct simulation", "l");
    legend->Draw();
        
    c1->SaveAs(Form("AngAcceptance.pdf"));
    delete c1;
    
    file1->Close();
    file2->Close();
}