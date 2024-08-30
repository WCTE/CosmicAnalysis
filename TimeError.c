R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)
//max_r and max_z is taken to be the inner boundary of the tank
#define max_r 159.627
#define max_z 141.5
#include "function.h"   //contains functions for main(muon) tracks selection and photon tracks reconstruction

void TimeError(const char* filename="./WCSimfitQunHitPreprocess.root") {
    //Open the file
    TFile* File = new TFile(filename, "read");
    if (!File->IsOpen()){
        std::cout << "Error, could not open the preprocess file: " << filename << "\n";
        return;
    }

    //Refer to code description of WCSimfitQunHitPreprocess.c for details of each branch
    TTree* EventTree = (TTree*)File->Get("EventTree");
    int THNum, DHNum;
    double fitQunVerTime, fitQunPos[3], fitQunDir[3], fitQunEntrance[3], fitQunExit[3], WCSimEntrance[3], WCSimStartPos[3], WCSimDir[3];
    EventTree->SetBranchAddress("THNum", &THNum);
    EventTree->SetBranchAddress("DHNum", &DHNum);
    EventTree->SetBranchAddress("fitQunVerTime", &fitQunVerTime);
    EventTree->SetBranchAddress("fitQunPos", fitQunPos);
    EventTree->SetBranchAddress("fitQunDir", fitQunDir);
    EventTree->SetBranchAddress("fitQunEntrance", fitQunEntrance);
    EventTree->SetBranchAddress("fitQunExit", fitQunExit);
    EventTree->SetBranchAddress("WCSimEntrance", WCSimEntrance);
    EventTree->SetBranchAddress("WCSimStartPos", WCSimStartPos);
    EventTree->SetBranchAddress("WCSimDir", WCSimDir);

    TTree* THTree = (TTree*)File->Get("THTree");
    int THPMT;
    double THTime;
    THTree->SetBranchAddress("THPMT", &THPMT);
    THTree->SetBranchAddress("THTime", &THTime);

    TTree* DHTree = (TTree*)File->Get("DHTree");
    int DHPMT;
    double DHTime;
    DHTree->SetBranchAddress("DHPMT", &DHPMT);
    DHTree->SetBranchAddress("DHTime", &DHTime);

    TTree* GeoTree = (TTree*)File->Get("GeoTree");
    double PMTOri[2014][3], PMTPos[2014][3];
    GeoTree->SetBranchAddress("PMTOri", PMTOri);
    GeoTree->SetBranchAddress("PMTPos", PMTPos);
    GeoTree->GetEntry(0);

    //Set up histogram
    TH1D* TrueTimeError = new TH1D("TrueTimeError","True time error;WCSim reconsturction - true hit time [ns];Number of event",1000,-10,10);
    TH1D* fitQunTimeError = new TH1D("fitQunTimeError","fitQun time error;fitQun reconsturction - digitized hit time [ns];Number of event",1000,-40,40);

    int TotalEventNum = EventTree->GetEntries();
    std::cout << "Number of events: " << TotalEventNum <<
        "\nNumber of true hits in all events: " << THTree->GetEntries() << 
        "\nNumber of digit hits in all events: " << DHTree->GetEntries() << "\n";

    int THTreeStartPos = 0, THTreeEndPos, DHTreeStartPos = 0, DHTreeEndPos;     //Starting and ending position of each event in the DHTree and THTree
    int TimeConstructionFailNum = 0;     //Number of main tracks that fail to reconstruct the entrance time (calculation gives a negative entrance time)
    int ignored_track = 0;      //Number of tracks filtered out by the selection() function
    
    //looping the event tree
    for (int i=0; i<TotalEventNum; i++) {
        EventTree->GetEntry(i);
        
        //selection of the fitQun(recontructed) main track. Note that the selection will only be apply on fitQun only but not WCSim(true)
        if (selection(fitQunEntrance, fitQunExit)) {
            THTreeStartPos += THNum;    //renew starting position of the THTree
            DHTreeStartPos += DHNum;    //renew starting position of the DHTree
            ignored_track++;
            continue;
        }

        double fitQunEntranceTime;    //Reconstructed time that the main track enters the tank
        //vertex inside the tank
        if ((fitQunPos[2]<max_z)&&((fitQunPos[0]*fitQunPos[0]+fitQunPos[1]*fitQunPos[1])<(max_r*max_r))) {
            fitQunEntranceTime = fitQunVerTime -
                (1e9)*sqrt((fitQunEntrance[0]-fitQunPos[0])*(fitQunEntrance[0]-fitQunPos[0]) +
                (fitQunEntrance[1]-fitQunPos[1])*(fitQunEntrance[1]-fitQunPos[1]) +
                (fitQunEntrance[2]-fitQunPos[2])*(fitQunEntrance[2]-fitQunPos[2]))/TMath::Ccgs();
        }
        //vertex outside the tank
        else {
            fitQunEntranceTime = fitQunVerTime +
                (1e9)*sqrt((fitQunEntrance[0]-fitQunPos[0])*(fitQunEntrance[0]-fitQunPos[0]) +
                (fitQunEntrance[1]-fitQunPos[1])*(fitQunEntrance[1]-fitQunPos[1]) +
                (fitQunEntrance[2]-fitQunPos[2])*(fitQunEntrance[2]-fitQunPos[2]))/TMath::Ccgs();
        }
        //reject negative entrance time
        if (fitQunEntranceTime<0) {
            THTreeStartPos += THNum;    //renew starting position of THTree
            DHTreeStartPos += DHNum;    //renew starting position of DHTree
            TimeConstructionFailNum++;
            continue;
        }

        const std::vector<std::vector<double>> fitQunPhoTrack = TrackReconstruction(fitQunEntrance, fitQunDir, PMTPos, PMTOri);     //2D array contains three parameters of a reconstruct photon track(from the fitQun vertex) for all 2014 PMTs. Check the code description of "function.h" for more details
        DHTreeEndPos = DHTreeStartPos + DHNum;      //renew ending position of DHTree
        //looping the THTree
        for (int j=DHTreeStartPos; j<DHTreeEndPos; j++) {
            DHTree->GetEntry(j);
            if (fitQunPhoTrack[DHPMT][0]<0) {continue;}     //Skipping PMTs which are not expected to be hit by Cherenkov photon
            fitQunTimeError->Fill(fitQunEntranceTime + fitQunPhoTrack[DHPMT][0] - DHTime);  //Fill time error histogram
        }
        DHTreeStartPos = DHTreeEndPos;      //renew starting position of DHTree

        double WCSimEntranceTime = (1e9)*sqrt((WCSimEntrance[0]-WCSimStartPos[0])*(WCSimEntrance[0]-WCSimStartPos[0]) +
            (WCSimEntrance[1]-WCSimStartPos[1])*(WCSimEntrance[1]-WCSimStartPos[1]) +
            (WCSimEntrance[2]-WCSimStartPos[2])*(WCSimEntrance[2]-WCSimStartPos[2]))/TMath::Ccgs();         //True time that the main track enters the tank

        const std::vector<std::vector<double>> TruePhoTrack = TrackReconstruction(WCSimEntrance, WCSimDir, PMTPos, PMTOri);     //2D array contains three parameters of a true photon track(from the WCSim vertex) for all PMTs
        THTreeEndPos = THTreeStartPos + THNum;      //renew ending position of THTree
        for (int j=THTreeStartPos; j<THTreeEndPos; j++) {
            THTree->GetEntry(j);
            if (TruePhoTrack[THPMT][0]<0) {continue;}     //Skipping PMTs which are not expected to be hit by Cherenkov photon
            TrueTimeError->Fill(WCSimEntranceTime + TruePhoTrack[THPMT][0] - THTime);  //Fill time error histogram
        }
        THTreeStartPos = THTreeEndPos;      //renew starting position of THTree
    }

    std::cout << "Number of tracks ignored: " << ignored_track <<
        "\nNumber of digitized hits that fail to reconstructed main track entering time: " << TimeConstructionFailNum << "\n";

    //Draw histogram
    TCanvas* c1 = new TCanvas();
    TrueTimeError->Draw();
    c1->SaveAs(Form("TrueTimeError.pdf"));
    delete TrueTimeError;
    fitQunTimeError->Draw();
    c1->SaveAs(Form("fitQunTimeError.pdf"));
    delete fitQunTimeError;

    File->Close();
    delete c1;
}