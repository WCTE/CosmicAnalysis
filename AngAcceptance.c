R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)
//max_r and max_z is taken to be the inner boundary of the tank
#define max_r 159.627
#define max_z 141.5
#include "function.h"

void AngAcceptance(const char* filename="./HitPreprocess.root") {
    TFile* File = new TFile(filename, "read");
    if (!File->IsOpen()){
        std::cout << "Error, could not open the preprocess file: " << filename << "\n";
        return;
    }
    TTree* EventTree = (TTree*)File->Get("EventTree");
    int DHNum, THNum;
    double fitQunVerTime, fitQunPos[3], fitQunDir[3], fitQunEntrance[3], fitQunExit[3], WCSimEntrance[3], WCSimStartPos[3], WCSimDir[3];
    EventTree->SetBranchAddress("DHNum", &DHNum);
    EventTree->SetBranchAddress("THNum", &THNum);
    EventTree->SetBranchAddress("fitQunVerTime", &fitQunVerTime);
    EventTree->SetBranchAddress("fitQunPos", fitQunPos);
    EventTree->SetBranchAddress("fitQunDir", fitQunDir);
    EventTree->SetBranchAddress("fitQunEntrance", fitQunEntrance);
    EventTree->SetBranchAddress("fitQunExit", fitQunExit);
    EventTree->SetBranchAddress("WCSimEntrance", WCSimEntrance);
    EventTree->SetBranchAddress("WCSimStartPos", WCSimStartPos);
    EventTree->SetBranchAddress("WCSimDir", WCSimDir);

    TTree* DHTree = (TTree*)File->Get("DHTree");
    int DHPMT;
    double DHTime, Charge;
    DHTree->SetBranchAddress("DHPMT", &DHPMT);
    DHTree->SetBranchAddress("DHTime", &DHTime);
    DHTree->SetBranchAddress("Charge", &Charge);

    TTree* THTree = (TTree*)File->Get("THTree");
    int THPMT, DHID;
    double THTime;
    THTree->SetBranchAddress("DHID", &DHID);
    THTree->SetBranchAddress("THPMT", &THPMT);
    THTree->SetBranchAddress("THTime", &THTime);

    TTree* GeoTree = (TTree*)File->Get("GeoTree");
    double PMTOri[2014][3], PMTPos[2014][3];
    GeoTree->SetBranchAddress("PMTOri", PMTOri);
    GeoTree->SetBranchAddress("PMTPos", PMTPos);
    GeoTree->GetEntry(0);

    //Set up histogram
    TH1D* fitQunDH = new TH1D("fitQunDH","Angular acceptance against cosine of the reconstructed incident angle;cosine;angular acceptance",10,0,1);
    TH1D* fitQunInAngle = new TH1D("fitQunInAngle","Cosine of the fitQun incident angle;cosine;Number of hit",10,0,1);
    TH1D* WCSimTH = new TH1D("WCSimTH","Angular acceptance against cosine of true incident angle;cosine;angular acceptance",10,0,1);
    TH1D* TrueInAngle = new TH1D("TrueInAngle","Cosine of the true incident angle;cosine;Number of true hits",10,0,1);
    TH1D* fitQunTH = new TH1D("fitQunTH","Angular acceptance against cosine of true incident angle;cosine;angular acceptance",10,0,1);

    TH1D* TrueBinHist[10];
    TH1D* fitQunBinHist[10];
    for (int i=0; i<10; i++) {
        TrueBinHist[i] = new TH1D(("True"+to_string(i)).c_str()," ",450,0,9000);
        fitQunBinHist[i] = new TH1D(("fitQun"+to_string(i)).c_str()," ",450,0,9000);
    }

    int TotalEventNum = EventTree->GetEntries();
    std::cout << "Number of events: " << TotalEventNum <<
        "\nNumber of digit hits in all events: " << DHTree->GetEntries() << "\n";

    int DHTreeStartPos = 0, DHTreeEndPos, THTreeStartPos = 0, THTreeEndPos;
    int TimeConstructionFailNum = 0, ignored_track = 0;
    for (int i=0; i<TotalEventNum; i++) {
        EventTree->GetEntry(i);
        std::cout << "\r" << i << std::flush;  //Enable this line to show the progress in runtime
        if (selection(fitQunEntrance, fitQunExit)) {
            DHTreeStartPos += DHNum;
            THTreeStartPos += THNum;
            ignored_track++;
            continue;
        }

        double fitQunEntranceTime;
        if ((fitQunPos[2]<max_z)&&((fitQunPos[0]*fitQunPos[0]+fitQunPos[1]*fitQunPos[1])<(max_r*max_r))) {
            fitQunEntranceTime = fitQunVerTime -
                (1e9)*sqrt((fitQunEntrance[0]-fitQunPos[0])*(fitQunEntrance[0]-fitQunPos[0]) +
                (fitQunEntrance[1]-fitQunPos[1])*(fitQunEntrance[1]-fitQunPos[1]) +
                (fitQunEntrance[2]-fitQunPos[2])*(fitQunEntrance[2]-fitQunPos[2]))/TMath::Ccgs();
        }
        else {
            fitQunEntranceTime = fitQunVerTime +
                (1e9)*sqrt((fitQunEntrance[0]-fitQunPos[0])*(fitQunEntrance[0]-fitQunPos[0]) +
                (fitQunEntrance[1]-fitQunPos[1])*(fitQunEntrance[1]-fitQunPos[1]) +
                (fitQunEntrance[2]-fitQunPos[2])*(fitQunEntrance[2]-fitQunPos[2]))/TMath::Ccgs();
        }
        if (fitQunEntranceTime<0) {
            DHTreeStartPos += DHNum;
            THTreeStartPos += THNum;
            TimeConstructionFailNum++;
            continue;
        }

        const std::vector<std::vector<double>> fitQunPhoTrack = TrackReconstruction(fitQunEntrance, fitQunDir, PMTPos, PMTOri);
        double PMTQ1[2014];
        for (int j=0; j<2014; j++) {
            if (fitQunPhoTrack[j][0]==-1) {continue;}
            PMTQ1[j] = 0;
            fitQunInAngle->Fill(fitQunPhoTrack[j][1]);
        }
        DHTreeEndPos = DHTreeStartPos + DHNum;
        for (int j=DHTreeStartPos; j<DHTreeEndPos; j++) {
            DHTree->GetEntry(j);
            if (fitQunPhoTrack[DHPMT][0]==-1) {continue;}
            double fitQunTimeError = fitQunEntranceTime + fitQunPhoTrack[DHPMT][0] - DHTime;
            if (fitQunTimeError<-2.1492) {continue;}
            PMTQ1[DHPMT] += Charge;
        }
        DHTreeStartPos = DHTreeEndPos;
        for (int j=0; j<2014; j++) {
            if (fitQunPhoTrack[j][0]==-1) {continue;}
            fitQunBinHist[(int)std::floor(fitQunPhoTrack[j][1]*10)]->Fill(PMTQ1[j]*fitQunPhoTrack[j][2]);
            if (PMTQ1[j]==0) {continue;}
            fitQunDH->Fill(fitQunPhoTrack[j][1],PMTQ1[j]*fitQunPhoTrack[j][2]);
        }
        
        double WCSimEntranceTime = (1e9)*sqrt((WCSimEntrance[0]-WCSimStartPos[0])*(WCSimEntrance[0]-WCSimStartPos[0]) +
            (WCSimEntrance[1]-WCSimStartPos[1])*(WCSimEntrance[1]-WCSimStartPos[1]) +
            (WCSimEntrance[2]-WCSimStartPos[2])*(WCSimEntrance[2]-WCSimStartPos[2]))/TMath::Ccgs();

        const std::vector<std::vector<double>> TruePhoTrack = TrackReconstruction(WCSimEntrance, WCSimDir, PMTPos, PMTOri);
        double PMTQ2[2014];
        for (int j=0; j<2014; j++) {
            if (TruePhoTrack[j][0]==-1) {continue;}
            PMTQ1[j] = 0;
            PMTQ2[j] = 0;
            TrueInAngle->Fill(TruePhoTrack[j][1]);
        }
        THTreeEndPos = THTreeStartPos + THNum;
        for (int j=THTreeStartPos; j<THTreeEndPos; j++) {
            THTree->GetEntry(j);
            double TrueTimeError = WCSimEntranceTime + TruePhoTrack[THPMT][0] - THTime;
            if (TrueTimeError<-1.7599) {continue;}
            if (TruePhoTrack[THPMT][0]!=(-1)) {PMTQ1[THPMT]++;}
            if (fitQunPhoTrack[THPMT][0]!=(-1)) {PMTQ2[THPMT]++;}
        }
        THTreeStartPos = THTreeEndPos;
        for (int j=0; j<2014; j++) {
            if (TruePhoTrack[j][0]==-1) {continue;}
            TrueBinHist[(int)std::floor(TruePhoTrack[j][1]*10)]->Fill(PMTQ1[j]*TruePhoTrack[j][2]);
            if (PMTQ1[j]!=0) {WCSimTH->Fill(TruePhoTrack[j][1],PMTQ1[j]*TruePhoTrack[j][2]);}
            if (PMTQ2[j]!=0) {fitQunTH->Fill(TruePhoTrack[j][1],PMTQ2[j]*TruePhoTrack[j][2]);}
        }
    }

    std::cout << "Number of tracks ignored: " << ignored_track << "\n" << 
        "Number of fitQun recontructed tracks which fails to calculate the entrance time: " << TimeConstructionFailNum << "\n";

    TCanvas* c1 = new TCanvas();
    TrueInAngle->Draw();
    c1->SaveAs(Form("TrueInAngle.pdf"));
    fitQunInAngle->Draw();
    c1->SaveAs(Form("fitQunInAngle.pdf"));

    int nbins = WCSimTH->GetNbinsX();
    for (int i=1; i<=nbins; i++) {
        double CurrentValue = WCSimTH->GetBinContent(i);
        WCSimTH->SetBinContent(i, CurrentValue/(TrueInAngle->GetBinContent(i)));
    }
    double NormalizationFactor = WCSimTH->GetBinContent(nbins);
    for (int i=1; i<=nbins; i++) {
        double CurrentValue = WCSimTH->GetBinContent(i);
        WCSimTH->SetBinContent(i, CurrentValue/NormalizationFactor);
    }

    for (int i=1; i<=nbins; i++) {
        double CurrentValue = fitQunDH->GetBinContent(i);
        fitQunDH->SetBinContent(i, CurrentValue/(fitQunInAngle->GetBinContent(i)));
    }
    NormalizationFactor = fitQunDH->GetBinContent(nbins);
    for (int i=1; i<=nbins; i++) {
        double CurrentValue = fitQunDH->GetBinContent(i);
        fitQunDH->SetBinContent(i, CurrentValue/NormalizationFactor);
    }

    for (int i=1; i<=nbins; i++) {
        double CurrentValue = fitQunTH->GetBinContent(i);
        fitQunTH->SetBinContent(i, CurrentValue/(fitQunInAngle->GetBinContent(i)));
    }
    NormalizationFactor = fitQunTH->GetBinContent(nbins);
    for (int i=1; i<=nbins; i++) {
        double CurrentValue = fitQunTH->GetBinContent(i);
        fitQunTH->SetBinContent(i, CurrentValue/NormalizationFactor);
    }

    for (int i=0; i<10; i++) {
        TrueBinHist[i]->SetStats(0);
        TrueBinHist[i]->SetLineColor(kRed);
        fitQunBinHist[i]->SetStats(0);
        fitQunBinHist[i]->SetLineColor(kBlue);
        std::string title = "Charge*distance distribution (Cosine of incident angle: [0." + std::to_string(i) + "-0." + std::to_string(i+1) + "])";
        TrueBinHist[i]->SetTitle(title.c_str());
        TrueBinHist[i]->SetXTitle("Charge * distance");
        TrueBinHist[i]->SetYTitle("Number of PMT hits");
        TrueBinHist[i]->Draw();
        fitQunBinHist[i]->Draw("same");
        TLegend* legend = new TLegend(0.75,0.75,0.9,0.9);
        legend->AddEntry(TrueBinHist[i], "True", "l");
        legend->AddEntry(fitQunBinHist[i], "fitQun", "l");
        legend->Draw();
        title = "BinHist[0." + std::to_string(i) + "-0." + std::to_string(i+1) + "].pdf";
        c1->SaveAs(title.c_str());
    }

    TFile* Hist = new TFile("AngAcceptance.root", "RECREATE");
    WCSimTH->Write();
    fitQunDH->Write();
    fitQunTH->Write();
    Hist->Close();

    delete c1;
    File->Close();
}
