//The code description "//g" below means that is the final value to be filled in to the tree

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)
//max_r and max_z is taken to be the inner boundary of the tank
#define max_r 159.627
#define max_z 141.5

bool extrapolation(const double fitQunPos[3], const double fitQunDir[3], double entrance[3], double exit[3], const double R=max_r, const double Z=max_z) {
    double z = fitQunPos[2];
    if (z<(-Z)) {return true;}      //select vertices with z coordinate not lower than the bottom of the tank
    double z_dir = fitQunDir[2];
    if (z_dir>0) {return true;}     //select vertices with downward travel direction
    double x = fitQunPos[0], y = fitQunPos[1];
    double x_dir = fitQunDir[0], y_dir = fitQunDir[1];
    double r = sqrt(x*x+y*y);
    double r_dir = sqrt(x_dir*x_dir+y_dir*y_dir);
    double vt;

    //check whether the vertex position is located in the tank
    if ((z<Z)&&(r<R)) {
        //calculate the exit point
        //assume the track exits through the bottom
        vt = (z-(-Z))/(-z_dir);
        exit[0] = x+x_dir*vt;
        exit[1] = y+y_dir*vt;
        exit[2] = -Z;
        //check if the track exits through the barrel
        if ((exit[0]*exit[0]+exit[1]*exit[1])>(R*R)) {
            double cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
            vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
            exit[0] = x+x_dir*vt;
            exit[1] = y+y_dir*vt;
            exit[2] = z+z_dir*vt;
        }

        //search for entrance point by reversing the direction
        x_dir = -x_dir;
        y_dir = -y_dir;
        z_dir = -z_dir;
        //calculate the entrance point
        //assume the track enters from the top
        vt = (Z-z)/z_dir;
        entrance[0] = x+x_dir*vt;
        entrance[1] = y+y_dir*vt;
        entrance[2] = Z;
        //check if the track enters from the barrel
        if ((entrance[0]*entrance[0]+entrance[1]*entrance[1])>(R*R)) {
            double cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
            vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
            entrance[0] = x+x_dir*vt;
            entrance[1] = y+y_dir*vt;
            entrance[2] = z+z_dir*vt;
        }
        return false;
    }

    //for vertices with radius larger than the tank radius
    if (r>R) {
        double cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);    //angle on the x-y plane between the travel direction and vertex position vectors
        if (acos(cos_theta)>asin(R/r)) {return true;}   //check whether the track enters the volume of the vertical extention of the tank
        //calculate entrance point
        //assume the track enters through the barrel
        vt = (r/r_dir)*(cos_theta-sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
        entrance[0] = x+x_dir*vt;
        entrance[1] = y+y_dir*vt;
        entrance[2] = z+z_dir*vt;
        if (entrance[2]<(-Z)) {return true;}    //check if the track miss the tank
        if (entrance[2]>Z) {
            //assume the track enters through the top
            vt = (z-Z)/(-z_dir);
            entrance[0] = x+x_dir*vt;
            entrance[1] = y+y_dir*vt;
            entrance[2] = Z;
            if ((entrance[0]*entrance[0]+entrance[1]*entrance[1])>(R*R)) {return true;}     //check if track enters through the top
        }

        //calculate exit position
        //assume the track exits through the bottom
        vt = (z-(-Z))/(-z_dir);
        exit[0] = x+x_dir*vt;
        exit[1] = y+y_dir*vt;
        exit[2] = -Z;
        //check if the track exits through the barrel
        if ((exit[0]*exit[0]+exit[1]*exit[1])>(R*R)) {
            vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
            exit[0] = x+x_dir*vt;
            exit[1] = y+y_dir*vt;
            exit[2] = z+z_dir*vt;
        }
        return false;
    }

    //for vertices on top of the tank with the vertex radius smaller than the tank radius
    if (z>Z) {
        //calculate the entrance point
        //assume the track enters from the top
        vt = (z-Z)/(-z_dir);
        entrance[0] = x+x_dir*vt;
        entrance[1] = y+y_dir*vt;
        entrance[2] = Z;
        if ((entrance[0]*entrance[0]+entrance[1]*entrance[1])>(R*R)) {return true;}     //check if the track enters from the top
        //calculate the exit point
        //assume the track exits through the bottom
        vt = (z-(-Z))/(-z_dir);
        exit[0] = x+x_dir*vt;
        exit[1] = y+y_dir*vt;
        exit[2] = -Z;
        //check if the particle exits through the barrel
        if ((exit[0]*exit[0]+exit[1]*exit[1])>(R*R)) {
            double cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
            vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
            exit[0] = x+x_dir*vt;
            exit[1] = y+y_dir*vt;
            exit[2] = z+z_dir*vt;
        }
        return false;
    }
    return true;
}

//selection of the fitQun vertex position
bool selection(const double VerPos[3], const double CutOff = max_r) {
    if ((VerPos[0]*VerPos[0]+VerPos[1]*VerPos[1])>=(CutOff*CutOff)) {return true;}
    return false;
}

//selection of the extrapolated entrance and exit position
bool selection(const double Entrance[3], const double Exit[3], const double CutOff = max_r) {
    if ((Entrance[0]*Entrance[0]+Entrance[1]*Entrance[1])>=(CutOff*CutOff)) {return true;}
    if ((Exit[0]*Exit[0]+Exit[1]*Exit[1])>=(CutOff*CutOff)) {return true;}
    return false;
}

void WCSimfitQunHitPreprocess(const char* fname="/work/kmtsui/wcte/cosmic/new_photocathode_model/fq_mu-_*.root", const char* wname="/work/kmtsui/wcte/cosmic/hk_flux/wcsim_mu-_*.root") {
    //Get the WCSim and fitQun files
    TChain* wt = new TChain("wcsimT");
    wt->Add(wname);
    TChain* ft = new TChain("fiTQun");
    ft->Add(fname);

    int nEntries = wt->GetEntries();
    std::cout << "Number of entries: " << nEntries << "\n";

    //Setup branches to read the fitQun and WCSim root files
    float fq1rpos[2][7][3];
    float fq1rdir[2][7][3];
    float fq1rt0[2][7];
    ft->SetBranchAddress("fq1rpos",fq1rpos);
    ft->SetBranchAddress("fq1rdir",fq1rdir);
    ft->SetBranchAddress("fq1rt0",fq1rt0);
    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    wt->SetBranchAddress("wcsimrootevent", &wcsimrootsuperevent);
    WCSimRootTrack* wcsimroottrack = new WCSimRootTrack();

    //Create the file "WCSimPreprocess.root"
    TFile* File = new TFile("WCSimfitQunHitPreprocess.root", "RECREATE");

    //Create the tree "EventTree"
    TTree* EventTree = new TTree("EventTree", "MainTrackData");
    int THNum, DHNum;
    double fitQunVerTime, fitQunPos[3], fitQunDir[3], fitQunEntrance[3], fitQunExit[3];
    double WCSimEntrance[3], WCSimExit[3], WCSimDir[3], WCSimStartPos[3];
    EventTree->Branch("THNum", &THNum);     //Number of true PMT hits in this event
    EventTree->Branch("DHNum", &DHNum);     //Number of digitized PMT hits in this event
    EventTree->Branch("fitQunVerTime", &fitQunVerTime);     //Reconstructed vertex time, Note that t=0 refers the time when the main track starts in the simulation
    EventTree->Branch("fitQunPos", fitQunPos, "fitQunPos[3]/D");     //Reconstructed vertex position
    EventTree->Branch("fitQunDir", fitQunDir, "fitQunDir[3]/D");     //Reconstructed main track direction
    EventTree->Branch("fitQunEntrance", fitQunEntrance, "fitQunEntrance[3]/D");     //Entrance point by extrapolating the reconstructed vertex
    EventTree->Branch("fitQunExit", fitQunExit, "fitQunExit[3]/D");     //Exit point by extrapolating the reconstructed vertex
    EventTree->Branch("WCSimEntrance", WCSimEntrance, "WCSimEntrance[3]/D");     //Entrance point of the true main track
    EventTree->Branch("WCSimExit", WCSimExit, "WCSimExit[3]/D");     //Exit point of the true main track
    EventTree->Branch("WCSimDir", WCSimDir, "WCSimDir[3]/D");     //Travel direction of the true main track
    EventTree->Branch("WCSimStartPos", WCSimStartPos, "WCSimStartPos[3]/D");     //Starting position of the main track particle in simulation

    //Create tree "GeoTree", which contains two 2D array carrying the position and orientation of all PMT. Note that the PMT with id=i will be saved in the i-1 th entry of the two arrays
    TTree* GeoTree = new TTree("GeoTree", "PMTGeometryData");
    double PMTOri[2014][3], PMTPos[2014][3];    //There are 2014 PMT in the detector
    GeoTree->Branch("PMTOri", PMTOri, "PMTOri[2014][3]/D");     //Array carrying the orientation of the PMTs
    GeoTree->Branch("PMTPos", PMTPos, "PMTPos[2014][3]/D");     //Array carrying the position of the PMTs
    //Get tree for geometry
    TChain* GeoT = new TChain("wcsimGeoT");
    GeoT->Add(wname);
    WCSimRootGeom* Geo = new WCSimRootGeom();
    GeoT->SetBranchAddress("wcsimrootgeom", &Geo);
    GeoT->GetEntry(0);
    if (Geo->GetWCNumPMT()!=2014) {
        std::cout << "Incorrect PMT numbers\n";
        return;
    }
    //looping for the PMTs
    for (int i=0;i<2014;i++) {
        PMTOri[i][0] = -(Geo->GetPMT(i)).GetOrientation(0);     //g
        PMTOri[i][1] = (Geo->GetPMT(i)).GetOrientation(2);     //g
        PMTOri[i][2] = (Geo->GetPMT(i)).GetOrientation(1);     //g
        PMTPos[i][0] = -(Geo->GetPMT(i)).GetPosition(0);     //g
        PMTPos[i][1] = (Geo->GetPMT(i)).GetPosition(2);     //g
        PMTPos[i][2] = (Geo->GetPMT(i)).GetPosition(1);     //g
    }
    GeoTree->Fill();    //fill the GeoTree with the new entry
    GeoTree->Write("", TObject::kOverwrite);  //Save the GeoTree in the file

    //Create true hit tree "THTree", each entry contains information of a true PMT hit
    TTree* THTree = new TTree("THTree", "PMTTrueHitData");
    int THPMT, DHID;
    double THTime, PhoStartPos[3], PhoEndPos[3], PhoDir[3], PhoStartDir[3], PhoEndDir[3];
    THTree->Branch("THPMT", &THPMT);    //Position of the hit PMT in the PMTOri and PMTPos array
    THTree->Branch("DHID", &DHID);      //Position of the corresponding digitized PMT hit of this true hit in DHTree 
    THTree->Branch("THTime", &THTime);      //True hit time
    THTree->Branch("PhoStartPos", PhoStartPos, "PhoStartPos[3]/D");     //Starting position of the photon
    THTree->Branch("PhoEndPos", PhoEndPos, "PhoEndPos[3]/D");     //Ending position of the photon
    THTree->Branch("PhoDir", PhoDir, "PhoDir[3]/D");     //Photon direction obtained by the starting and ending position
    THTree->Branch("PhoStartDir", PhoStartDir, "PhoStartDir[3]/D");     //Photon direction at the starting position
    THTree->Branch("PhoEndDir", PhoEndDir, "PhoEndDir[3]/D");     //Photon direction at the ending position

    //Create digitized hit tree "DHTree", contains the digitized hit information
    TTree* DHTree = new TTree("DHTree", "PMTDigitHitData");
    int DHPMT;
    double DHTime, Charge;
    DHTree->Branch("DHPMT", &DHPMT);    //Position of the hit PMT in the PMTOri and PMTPos array
    DHTree->Branch("DHTime", &DHTime);      //Digit hit time
    DHTree->Branch("Charge", &Charge);      //Charge (Reconstructed number of photon hit of a PMT)

    int ignored_track = 0;      //record the number of track being filtered out
    for (int i=0; i<nEntries; i++) {
        std::cout << "\r" << i << std::flush;  //Enable this line to show the progress in runtime
        ft->GetEntry(i);
        fitQunPos[0] = -fq1rpos[0][2][0];   //g
        fitQunPos[1] = fq1rpos[0][2][2];   //g
        fitQunPos[2] = fq1rpos[0][2][1];   //g
        //selection on the main track by the fitQun vertex position
        if (selection(fitQunPos)) {
            ignored_track++;
            continue;
        }
        fitQunDir[0] = -fq1rdir[0][2][0];   //g
        fitQunDir[1] = fq1rdir[0][2][2];   //g
        fitQunDir[2] = fq1rdir[0][2][1];   //g

        //extrapolate the fitQun vertex to obtain the entrance and exit position (fitQunEntrance[3] and fitQunExit[3])
        if (extrapolation(fitQunPos, fitQunDir, fitQunEntrance, fitQunExit)) {
            ignored_track++;
            continue;
        }   //g for fitQunEntrance[3] and fitQunExit[3]

        //selection on the main track by the fitQun entrance and exit position
        if (selection(fitQunEntrance, fitQunExit)) {
            ignored_track++;
            continue;
        }
        fitQunVerTime = fq1rt0[0][2];   //g

        wt->GetEntry(i);
        WCSimRootTrigger* wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
        //search for main track
        int ntrack = wcsimrootevent->GetNtrack_slots();
        for (int itrack=0; itrack<ntrack; itrack++) {
            wcsimroottrack = dynamic_cast<WCSimRootTrack*> (wcsimrootevent->GetTracks()->At(itrack));
            if (wcsimroottrack->GetId()==1) {break;}
        }
        int nCrossing = wcsimroottrack->GetBoundaryPoints().size();
        if (nCrossing!=4) {
            ignored_track++;
            continue;
        } //ignoring tracks that did not enter or terminated in the tank

        //calculate WCSim entrance and exit points
        WCSimEntrance[0] = -wcsimroottrack->GetBoundaryPoints().at(1).at(0)/10;    //g
        WCSimEntrance[1] = wcsimroottrack->GetBoundaryPoints().at(1).at(2)/10;    //g
        WCSimEntrance[2] = wcsimroottrack->GetBoundaryPoints().at(1).at(1)/10;    //g
        WCSimExit[0] = -wcsimroottrack->GetBoundaryPoints().at(2).at(0)/10;    //g
        WCSimExit[1] = wcsimroottrack->GetBoundaryPoints().at(2).at(2)/10;    //g
        WCSimExit[2] = wcsimroottrack->GetBoundaryPoints().at(2).at(1)/10;    //g

        //calculate the WCSim travel direction
        WCSimDir[0] = WCSimExit[0] - WCSimEntrance[0];
        WCSimDir[1] = WCSimExit[1] - WCSimEntrance[1];
        WCSimDir[2] = WCSimExit[2] - WCSimEntrance[2];
        double WCSimDist = sqrt(WCSimDir[0]*WCSimDir[0]+WCSimDir[1]*WCSimDir[1]+WCSimDir[2]*WCSimDir[2]);
        WCSimDir[0] = WCSimDir[0]/WCSimDist;    //g
        WCSimDir[1] = WCSimDir[1]/WCSimDist;    //g
        WCSimDir[2] = WCSimDir[2]/WCSimDist;    //g
        WCSimStartPos[0] = -wcsimroottrack->GetStart(0);    //g
        WCSimStartPos[1] = wcsimroottrack->GetStart(2);    //g
        WCSimStartPos[2] = wcsimroottrack->GetStart(1);    //g
        THNum = wcsimrootevent->GetNcherenkovhittimes();

        TClonesArray* WRCDHArray = wcsimrootevent->GetCherenkovDigiHits();    //Array of WCSimRootCherenkovDigiHit pointer
        DHNum = WRCDHArray->GetEntriesFast();
        for (int j=0; j<DHNum; j++) {
            //ignore extremely long digit time (>200 ns)
            if (((WCSimRootCherenkovDigiHit*) WRCDHArray->At(j))->GetT()>200) {
                WRCDHArray->RemoveAt(j);
            }
        }
        WRCDHArray->Compress();
        DHNum = WRCDHArray->GetEntriesFast();    //g
        //loop the array and fill the DHTree
        for (int j=0; j<DHNum; j++) {
            WCSimRootCherenkovDigiHit* WRCDH = (WCSimRootCherenkovDigiHit*) WRCDHArray->At(j);
            DHPMT = WRCDH->GetTubeId()-1;   //g
            DHTime = WRCDH->GetT();     //g
            Charge = WRCDH->GetQ();     //g
            DHTree->Fill();
        }

        //loop for the PMTs with true hit in this event
        int THPMTNum = wcsimrootevent->GetNcherenkovhits();
        for (int j=0; j<THPMTNum; j++) {
            WCSimRootCherenkovHit* WRCH = (WCSimRootCherenkovHit*) wcsimrootevent->GetCherenkovHits()->At(j);
            THPMT = WRCH->GetTubeID()-1;    //g

            int FirstHitPos = WRCH->GetTotalPe(0), LastHitPos = WRCH->GetTotalPe(0)+WRCH->GetTotalPe(1)-1;     //The position on the Array of WCSimRootCherenkovHitTime of the first hit and the last hit on this PMT
            //loop for the true hits of a PMT in an event and fill the THTree
            for (int k=FirstHitPos; k<=LastHitPos; k++) {
                WCSimRootCherenkovHitTime* WRCHT = (WCSimRootCherenkovHitTime*) wcsimrootevent->GetCherenkovHitTimes()->At(k);
                if (WRCHT->GetParentSavedTrackID()==(-1)) {
                    THNum--;
                    continue;
                }   //ignore dark noise

                DHID = -1;  //The DHTree position will be set to -1 if the true hit does not generate a digit hit
                for (int n=0; n<DHNum; n++) {
                    WCSimRootCherenkovDigiHit* WRCDH = (WCSimRootCherenkovDigiHit*) WRCDHArray->At(n);
                    std::vector<int> fPhotonIds = WRCDH->GetPhotonIds();
                    for (int m : fPhotonIds) {
                        if (m!=k) {continue;}
                        DHID = DHTree->GetEntries()-1 - (DHNum-1-n);    //g
                        break;
                    }
                    if (DHID!=(-1)) {break;}
                }

                THTime = WRCHT->GetTruetime();    //g
                PhoStartPos[0] = -WRCHT->GetPhotonStartPos(0)/10;    //g
                PhoStartPos[1] = WRCHT->GetPhotonStartPos(2)/10;    //g
                PhoStartPos[2] = WRCHT->GetPhotonStartPos(1)/10;    //g
                PhoEndPos[0] = -WRCHT->GetPhotonEndPos(0)/10;    //g
                PhoEndPos[1] = WRCHT->GetPhotonEndPos(2)/10;    //g
                PhoEndPos[2] = WRCHT->GetPhotonEndPos(1)/10;    //g

                PhoStartDir[0] = -WRCHT->GetPhotonStartDir(0);    //g
                PhoStartDir[1] = WRCHT->GetPhotonStartDir(2);    //g
                PhoStartDir[2] = WRCHT->GetPhotonStartDir(1);    //g

                PhoEndDir[0] = -WRCHT->GetPhotonEndDir(0);    //g
                PhoEndDir[1] = WRCHT->GetPhotonEndDir(2);    //g
                PhoEndDir[2] = WRCHT->GetPhotonEndDir(1);    //g

                PhoDir[0] = PhoEndPos[0] - PhoStartPos[0];
                PhoDir[1] = PhoEndPos[1] - PhoStartPos[1];
                PhoDir[2] = PhoEndPos[2] - PhoStartPos[2];
                double PhoDist = sqrt(PhoDir[0]*PhoDir[0] + PhoDir[1]*PhoDir[1] + PhoDir[2]*PhoDir[2]);
                PhoDir[0] = PhoDir[0]/PhoDist;    //g
                PhoDir[1] = PhoDir[1]/PhoDist;    //g
                PhoDir[2] = PhoDir[2]/PhoDist;    //g

                THTree->Fill();    //fill the THTree with the new entry
            }
        }    //g for THNum
        EventTree->Fill();  //fill the EventTree with the new entry
    }

    std::cout << "Number of tracks ignored: " << ignored_track << "\n";
    EventTree->Write("", TObject::kOverwrite);  //Save the EventTree in the file
    DHTree->Write("", TObject::kOverwrite);  //Save the DHTree in the file
    THTree->Write("", TObject::kOverwrite);  //Save the THTree in the file
    File->Close();
}