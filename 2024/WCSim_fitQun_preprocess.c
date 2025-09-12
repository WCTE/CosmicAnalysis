//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542
R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

//function for extrapolation of the reconstructed vertex, return true when the track formed by extrapolation of vertex does not across the tank
bool extrapolation(float fq1rpos[100][7][3], float fq1rdir[100][7][3], double entrance[3], double exit[3], double R=max_r, double Z=max_z) {
    double z = fq1rpos[0][2][1];
    if (z<(-Z)) {return true;}      //select vertices with z coordinate not lower than the bottom of the tank
    double z_dir = fq1rdir[0][2][1];
    if (z_dir>0) {return true;}     //select vertices with downward travel direction
    double x = -fq1rpos[0][2][0], y = fq1rpos[0][2][2];
    double x_dir = -fq1rdir[0][2][0], y_dir = fq1rdir[0][2][2];
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

void WCSim_fitQun_preprocess(const char* fname="/work/kmtsui/wcte/cosmic/new_photocathode_model/fq_mu-_*.root", const char* wname="/work/kmtsui/wcte/cosmic/new_photocathode_model/wcsim_mu-_*.root") {
    //Get the files
    TChain* ft = new TChain("fiTQun");
    ft->Add(fname);
    TChain* wt = new TChain("wcsimT");
    wt->Add(wname);

    int nEntries = ft->GetEntries();
    if (nEntries!=(wt->GetEntries())) {
        std::cout << "Number of events of WCSim and fitQun not equal\n";
        return;
    }
    else {std::cout << "Number of entries: " << nEntries << "\n";}

    //Set up branches for fitQun data
    float fq1rpos[100][7][3];
    float fq1rdir[100][7][3];
    float fq1rnll[100][7];
    ft->SetBranchAddress("fq1rpos",fq1rpos);
    ft->SetBranchAddress("fq1rdir",fq1rdir);
    ft->SetBranchAddress("fq1rnll",fq1rnll);

    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    wt->SetBranchAddress("wcsimrootevent", &wcsimrootsuperevent);
    WCSimRootTrack* wcsimroottrack = new WCSimRootTrack();

    //Create the file "WCSim_fitQun_preprocess.root"
    TFile* file = new TFile("WCSim_fitQun_preprocess.root", "RECREATE");
    TTree* tree = new TTree("WCSim_and_fitQun", "Preprocess_data");

    //Set up branch addresses for the data to be saved
    double WCSimEntrance[3], WCSimExit[3], WCSimDir[3], fitQunPos[3], fitQunEntrance[3], fitQunExit[3], fitQunDir[3];
    double WCSimTotalQ, WCSimDist, fitQunDist, fitQunLikelihood, eta;
    tree->Branch("WCSimEntrance", WCSimEntrance, Form("WCSimEntrance[%d]/D", 3));   //WCSim entrance point
    tree->Branch("WCSimExit", WCSimExit, Form("WCSimExit[%d]/D", 3));   //WCSim exit point
    tree->Branch("WCSimDir", WCSimDir, Form("WCSimDir[%d]/D", 3));  //WCSim travel direction
    tree->Branch("fitQunPos", fitQunPos, Form("fitQunPos[%d]/D", 3));   //fitQun vertex position
    tree->Branch("fitQunEntrance", fitQunEntrance, Form("fitQunEntrance[%d]/D", 3));    //fitQun extrapolated entrance point
    tree->Branch("fitQunExit", fitQunExit, Form("fitQunExit[%d]/D", 3));    //fitQun extrapolated exit point
    tree->Branch("fitQunDir", fitQunDir, Form("fitQunDir[%d]/D", 3));   //fitQun travel direction
    tree->Branch("WCSimTotalQ", &WCSimTotalQ);  //WCSim total PMT charge
    tree->Branch("WCSimDist", &WCSimDist);  //WCSim travel distance
    tree->Branch("fitQunDist", &fitQunDist);    //fitQun travel distance
    tree->Branch("fitQunLikelihood", &fitQunLikelihood);    //fitQun likelihood
    tree->Branch("eta", &eta);  //angle between the fitQun and WCSim travel direction

    for (int i=0; i<nEntries; i++) {
        //fitQun part
        ft->GetEntry(i);
        if (extrapolation(fq1rpos, fq1rdir, fitQunEntrance, fitQunExit)) {continue;}    //perform extrapolation

        //WCSim part
        delete wcsimrootsuperevent;
        wcsimrootsuperevent = 0;  // EXTREMELY IMPORTANT
        wt->GetEntry(i);
        WCSimRootTrigger* wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
        //search for main track
        int ntrack = wcsimrootevent->GetNtrack_slots();
        for (int itrack=0; itrack<ntrack; itrack++) {
            wcsimroottrack = dynamic_cast<WCSimRootTrack*> (wcsimrootevent->GetTracks()->At(itrack));
            if (wcsimroottrack->GetId()==1) {break;}
        }
        int nCrossing = wcsimroottrack->GetBoundaryPoints().size();
        if (nCrossing<4) {continue;}    //ignoring tracks that did not enter or terminated in the tank

        //calculate and save data
        //rotation
        fitQunPos[0] = -fq1rpos[0][2][0];
        fitQunPos[1] = fq1rpos[0][2][2];
        fitQunPos[2] = fq1rpos[0][2][1];
        fitQunDir[0] = -fq1rdir[0][2][0];
        fitQunDir[1] = fq1rdir[0][2][2];
        fitQunDir[2] = fq1rdir[0][2][1];
        fitQunDist = sqrt((fitQunEntrance[0]-fitQunExit[0])*(fitQunEntrance[0]-fitQunExit[0])+(fitQunEntrance[1]-fitQunExit[1])*(fitQunEntrance[1]-fitQunExit[1])+(fitQunEntrance[2]-fitQunExit[2])*(fitQunEntrance[2]-fitQunExit[2]));
        fitQunLikelihood = fq1rnll[0][2];

        //rotation
        WCSimEntrance[0] = -wcsimroottrack->GetBoundaryPoints().at(0).at(0)/10;
        WCSimEntrance[1] = wcsimroottrack->GetBoundaryPoints().at(0).at(2)/10;
        WCSimEntrance[2] = wcsimroottrack->GetBoundaryPoints().at(0).at(1)/10;
        WCSimExit[0] = -wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(0)/10;
        WCSimExit[1] = wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(2)/10;
        WCSimExit[2] = wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(1)/10;

        WCSimDir[0] = WCSimExit[0] - WCSimEntrance[0];
        WCSimDir[1] = WCSimExit[1] - WCSimEntrance[1];
        WCSimDir[2] = WCSimExit[2] - WCSimEntrance[2];
        WCSimDist = sqrt(WCSimDir[0]*WCSimDir[0]+WCSimDir[1]*WCSimDir[1]+WCSimDir[2]*WCSimDir[2]);
        WCSimDir[0] = WCSimDir[0]/WCSimDist;
        WCSimDir[1] = WCSimDir[1]/WCSimDist;
        WCSimDir[2] = WCSimDir[2]/WCSimDist;

        eta = acos(fitQunDir[0]*WCSimDir[0]+fitQunDir[1]*WCSimDir[1]+fitQunDir[2]*WCSimDir[2]);
        
        //Sum the charge of each PMTs and get the total
        WCSimTotalQ = 0;
        int nhits = wcsimrootevent->GetNcherenkovdigihits();
        for (int i=0; i<nhits; i++) {
            WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit = (WCSimRootCherenkovDigiHit*) (wcsimrootevent->GetCherenkovDigiHits())->At(i);
            WCSimTotalQ += wcsimrootcherenkovdigihit->GetQ();
        }

        tree->Fill();   //fill the tree with the new entry
    }

    tree->Write();  //Save the tree in the file
    file->Close();
    ft->Reset();
    wt->Reset();
}