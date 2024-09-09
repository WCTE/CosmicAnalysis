R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)
//max_r and max_z is taken to be the inner boundary of the tank
#define max_r 159.627
#define max_z 141.5
#define WaterCherenkovAng 0.719887812
#define WaterRefractiveIndex 1.33

//return true when the main track has entrance or exit radius larger than the cutoff, otherwise return false
//EntranceCut and ExitCut are the cutoff radius of the entrance and exit point
bool selection(const double Entrance[3], const double Exit[3], const double EntranceCut = 120, const double ExitCut = 145) {
    if ((Entrance[0]*Entrance[0]+Entrance[1]*Entrance[1])>(EntranceCut*EntranceCut)) {return true;}
    if ((Exit[0]*Exit[0]+Exit[0]*Exit[0])>(ExitCut*ExitCut)) {return true;}
    return false;
}

//the function takes in 4 arrays to reconstruct the track of the emitted Cherenkov photon from the main track to each PMT
//the 4 arguments are, the main track entrance position and direction, as well as all PMTs' position and orientation
//the function returns the hit time, incident angle and photon emission position as an array of that PMT. If reconstruction fails, it returns {-1}
std::vector<std::vector<double>> TrackReconstruction(const double VerPos[3], const double VerDir[3], const double PMTPos[2014][3], const double PMTOri[2014][3]) {
    std::vector<std::vector<double>> PMTReTrack;    //array to be returned
    double VerToPMTDir[3], PhoDir[3], VerToPMTDist, VerToPhoDist, PhoToPMTDist, theta, HitTime, CosInAngle;
    
    //looping all the PMTs
    for (int i=0; i<2014; i++) {
        //VerToPMTDir[3] is a unit vector pointing from the entrance position to the PMT
        VerToPMTDir[0] = PMTPos[i][0] - VerPos[0];
        VerToPMTDir[1] = PMTPos[i][1] - VerPos[1];
        VerToPMTDir[2] = PMTPos[i][2] - VerPos[2];
        VerToPMTDist = sqrt(VerToPMTDir[0]*VerToPMTDir[0] + VerToPMTDir[1]*VerToPMTDir[1] + VerToPMTDir[2]*VerToPMTDir[2]);     //distance from the entrance position of the main track to the target PMT
        VerToPMTDir[0] = VerToPMTDir[0]/VerToPMTDist;
        VerToPMTDir[1] = VerToPMTDir[1]/VerToPMTDist;
        VerToPMTDir[2] = VerToPMTDir[2]/VerToPMTDist;

        //theta is the angle between the entrance direction vector and the VerToPMTDir[3] vector
        theta = VerToPMTDir[0]*VerDir[0] + VerToPMTDir[1]*VerDir[1] + VerToPMTDir[2]*VerDir[2];     //This is the value of cos(theta)
        if (fabs(theta)>1) {theta = 1;}  //avoid numerial percision problem
        theta = acos(theta);
        //Reconstruction fail when theta > emission angle of the Cherenkov photon
        if (theta>=WaterCherenkovAng) {
            PMTReTrack.push_back({-1});
            continue;
        }
        VerToPhoDist = sin(WaterCherenkovAng - theta)*VerToPMTDist/sin(TMath::Pi()-WaterCherenkovAng);    //distance from the entrance position to the Cherenkov photon emission position
        //PhoDir[3] is the direction of the Cherenkov photon
        PhoDir[0] = PMTPos[i][0] - (VerPos[0]+VerDir[0]*VerToPhoDist);
        PhoDir[1] = PMTPos[i][1] - (VerPos[1]+VerDir[1]*VerToPhoDist);
        PhoDir[2] = PMTPos[i][2] - (VerPos[2]+VerDir[2]*VerToPhoDist);
        PhoToPMTDist = sqrt(PhoDir[0]*PhoDir[0] + PhoDir[1]*PhoDir[1] + PhoDir[2]*PhoDir[2]);    //Travel distance of the Cherenkov photon to the targeted PMT
        PhoDir[0] /= PhoToPMTDist;
        PhoDir[1] /= PhoToPMTDist;
        PhoDir[2] /= PhoToPMTDist;
        CosInAngle = -PMTOri[i][0]*PhoDir[0]-PMTOri[i][1]*PhoDir[1]-PMTOri[i][2]*PhoDir[2];     //Cosine of the incident angle
        //ignored the case that a photon enters a PMT at the back
        if (CosInAngle<0) {
            PMTReTrack.push_back({-1});
            continue;
        }
        HitTime = (1e9)*(VerToPhoDist + PhoToPMTDist*WaterRefractiveIndex)/TMath::Ccgs();   //Hit time of the photon, where t=0 is the time when the main track enters the tank
        //push back the array containing the information of the photon reconstruction of this PMT
        PMTReTrack.push_back({HitTime, CosInAngle, VerPos[0]+VerDir[0]*VerToPhoDist, VerPos[1]+VerDir[1]*VerToPhoDist, VerPos[2]+VerDir[2]*VerToPhoDist});
    }
    return PMTReTrack;
}

void AngPosError(const char* filename="./WCSimfitQunHitPreprocess.root") {
    TFile* File = new TFile(filename, "read");
    if (!File->IsOpen()){
        std::cout << "Error, could not open the preprocess file: " << filename << "\n";
        return;
    }

    //Refer to code description of WCSimfitQunHitPreprocess.c for details of each branch
    TTree* EventTree = (TTree*)File->Get("EventTree");
    int THNum;
    double fitQunVerTime, fitQunPos[3], fitQunDir[3], fitQunEntrance[3], fitQunExit[3], WCSimStartPos[3], WCSimEntrance[3], WCSimDir[3];
    EventTree->SetBranchAddress("THNum", &THNum);
    EventTree->SetBranchAddress("fitQunVerTime", &fitQunVerTime);
    EventTree->SetBranchAddress("fitQunPos", fitQunPos);
    EventTree->SetBranchAddress("fitQunDir", fitQunDir);
    EventTree->SetBranchAddress("fitQunEntrance", fitQunEntrance);
    EventTree->SetBranchAddress("fitQunExit", fitQunExit);
    EventTree->SetBranchAddress("WCSimStartPos", WCSimStartPos);
    EventTree->SetBranchAddress("WCSimEntrance", WCSimEntrance);
    EventTree->SetBranchAddress("WCSimDir", WCSimDir);

    TTree* THTree = (TTree*)File->Get("THTree");
    int THPMT, DHID;
    double THTime, PhoDir[3], PhoStartPos[3];
    THTree->SetBranchAddress("THPMT", &THPMT);
    THTree->SetBranchAddress("DHID", &DHID);
    THTree->SetBranchAddress("THTime", &THTime);
    THTree->SetBranchAddress("PhoDir", PhoDir);
    THTree->SetBranchAddress("PhoStartPos", PhoStartPos);

    TTree* DHTree = (TTree*)File->Get("DHTree");
    double DHTime;
    DHTree->SetBranchAddress("DHTime", &DHTime);

    TTree* GeoTree = (TTree*)File->Get("GeoTree");
    double PMTOri[2014][3], PMTPos[2014][3];
    GeoTree->SetBranchAddress("PMTOri", PMTOri);
    GeoTree->SetBranchAddress("PMTPos", PMTPos);
    GeoTree->GetEntry(0);

    //Set up histogram, details refer to Jimmy_README.md
    TH1D* fitQunAngError = new TH1D("fitQunAngError", "Error of fitQun reconstructed incident angle;cosine of fitQun incident angle - true incident angle;Number of true hits",100,-0.5,0.5);
    TH1D* fitQunPosError = new TH1D("fitQunPosError", "Error of fitQun reconstructed photon emission position;[cm];Number of true hits",500,0,100);
    TH1D* WCSimAngError = new TH1D("WCSimAngError", "Error of WCSim reconstructed incident angle;cosine of WCSim incident angle - true incident angle;Number of true hits",100,-0.5,0.5);
    TH1D* WCSimPosError = new TH1D("WCSimPosError", "Error of WCSim reconstructed photon emission position;[cm];Number of true hits",500,0,100);
    int TotalEventNum = EventTree->GetEntries();
    std::cout << "Number of events: " << TotalEventNum << "\n";

    int THTreeStartPos = 0, THTreeEndPos;   //Record the starting and ending position of true hits of a particular event in the THTree
    int ignored_track = 0;  //Number of main track ignored in selection

    //loop for each event
    for (int i=0; i<TotalEventNum; i++) {
        EventTree->GetEntry(i);
        std::cout << "\r" << i << std::flush;  //Enable this line to show the progress in runtime

        //select the main track base on the entrance and exit cutoff
        if (selection(fitQunEntrance, fitQunExit)) {
            THTreeStartPos += THNum;    //THTree starting position of the next event
            ignored_track++;
            continue;
        }

        //calculate the time when the main track enters the tank for later hit time calculation
        double fitQunEntranceTime;
        //for fitQun vertex inside the tank
        if ((fitQunPos[2]<max_z)&&((fitQunPos[0]*fitQunPos[0]+fitQunPos[1]*fitQunPos[1])<(max_r*max_r))) {
            fitQunEntranceTime = fitQunVerTime -
                (1e9)*sqrt((fitQunEntrance[0]-fitQunPos[0])*(fitQunEntrance[0]-fitQunPos[0]) +
                (fitQunEntrance[1]-fitQunPos[1])*(fitQunEntrance[1]-fitQunPos[1]) +
                (fitQunEntrance[2]-fitQunPos[2])*(fitQunEntrance[2]-fitQunPos[2]))/TMath::Ccgs();
        }
        //for fitQun vertex outside the tank
        else {
            fitQunEntranceTime = fitQunVerTime +
                (1e9)*sqrt((fitQunEntrance[0]-fitQunPos[0])*(fitQunEntrance[0]-fitQunPos[0]) +
                (fitQunEntrance[1]-fitQunPos[1])*(fitQunEntrance[1]-fitQunPos[1]) +
                (fitQunEntrance[2]-fitQunPos[2])*(fitQunEntrance[2]-fitQunPos[2]))/TMath::Ccgs();
        }

        //entrance time of main track according to WCSim
        double WCSimEntranceTime = (1e9)*sqrt((WCSimEntrance[0]-WCSimStartPos[0])*(WCSimEntrance[0]-WCSimStartPos[0]) +
            (WCSimEntrance[1]-WCSimStartPos[1])*(WCSimEntrance[1]-WCSimStartPos[1]) +
            (WCSimEntrance[2]-WCSimStartPos[2])*(WCSimEntrance[2]-WCSimStartPos[2]))/TMath::Ccgs();

        //recontruction of the Cherenkov photon track by the function TrackReconstruction()
        const std::vector<std::vector<double>> fitQunPhoTrack = TrackReconstruction(fitQunEntrance, fitQunDir, PMTPos, PMTOri);    //vector carrying information of the Cherenkov photon track recontructed from fitQun vertex
        const std::vector<std::vector<double>> TruePhoTrack = TrackReconstruction(WCSimEntrance, WCSimDir, PMTPos, PMTOri);    //vector carrying information of the Cherenkov photon track recontructed from WCSim

        //looping the true hits of this event
        THTreeEndPos = THTreeStartPos + THNum;
        for (int j=THTreeStartPos; j<THTreeEndPos; j++) {
            THTree->GetEntry(j);
            double TrueInAng = -PhoDir[0]*PMTOri[THPMT][0]-PhoDir[1]*PMTOri[THPMT][1]-PhoDir[2]*PMTOri[THPMT][2];   //Calculate the true incident angle
            if (TruePhoTrack[THPMT][0]!=(-1)) {
                double TrueTimeError = WCSimEntranceTime + TruePhoTrack[THPMT][0] - THTime;
                //timing cutoff
                if (TrueTimeError>-1.7589) {
                    WCSimAngError->Fill(TruePhoTrack[THPMT][1] - TrueInAng);    //WCSim reconstructed incident angle - true incident angle
                    WCSimPosError->Fill(sqrt((TruePhoTrack[THPMT][2]-PhoStartPos[0])*(TruePhoTrack[THPMT][2]-PhoStartPos[0]) +
                        (TruePhoTrack[THPMT][3]-PhoStartPos[1])*(TruePhoTrack[THPMT][3]-PhoStartPos[1]) +
                        (TruePhoTrack[THPMT][4]-PhoStartPos[2])*(TruePhoTrack[THPMT][4]-PhoStartPos[2])));  //3D distance between the reconstructed photon emission position and the true photon emission position
                }
            }
            //if the true hit generates a digit hit, calculate the angular and position error of the fitQun reconstructed photon track
            if (fitQunPhoTrack[THPMT][0]!=(-1) && DHID!=(-1)) {
                DHTree->GetEntry(DHID);     //get the corresponding digit hit of the true hit
                double fitQunTimeError = fitQunEntranceTime + fitQunPhoTrack[THPMT][0] - DHTime;
                //timing cutoff
                if (fitQunTimeError>-2.1753) {
                    fitQunAngError->Fill(fitQunPhoTrack[THPMT][1] - TrueInAng);    //fitQun reconstructed incident angle - true incident angle
                    fitQunPosError->Fill(sqrt((fitQunPhoTrack[THPMT][2]-PhoStartPos[0])*(fitQunPhoTrack[THPMT][2]-PhoStartPos[0]) +
                        (fitQunPhoTrack[THPMT][3]-PhoStartPos[1])*(fitQunPhoTrack[THPMT][3]-PhoStartPos[1]) +
                        (fitQunPhoTrack[THPMT][4]-PhoStartPos[2])*(fitQunPhoTrack[THPMT][4]-PhoStartPos[2])));
                }
            }
        }
        THTreeStartPos = THTreeEndPos;      //Renew the starting and ending position for the next event
    }

    std::cout << "Number of tracks ignored: " << ignored_track << "\n";
    
    TCanvas* c1 = new TCanvas();
    fitQunAngError->Draw();
    c1->SaveAs(Form("fitQunAngError.pdf"));
    fitQunPosError->Draw();
    c1->SaveAs(Form("fitQunPosError.pdf"));
    WCSimAngError->Draw();
    c1->SaveAs(Form("WCSimAngError.pdf"));
    WCSimPosError->Draw();
    c1->SaveAs(Form("WCSimPosError.pdf"));

    delete c1;
    File->Close();
}