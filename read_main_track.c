//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

//iEvent is the event(in the given file) to be read
void read_main_track(int iEvent=0, const char *filename="/work/kmtsui/wcte/cosmic/buggy/wcsim_mu-_0000.root") {

    //Get file
    TFile* file = new TFile(filename, "read");
    TTree* tree = (TTree*)file->Get("wcsimT");
    
    //Check whether the given event exist
    std::cout << "Target event: " << iEvent << "\n";
    int nEvent = tree->GetEntries();    //number of event contain in the file
    std::cout << "Total number of events: " << nEvent << "\n";
    if (iEvent >= nEvent) {
        std::cout << "Event " << iEvent << " does not exist\n";
        return;
    }

    //Construct the event branch
    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    tree->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent);
    
    //Fill the event branch with the event indicated by iEvent
    tree->GetEntry(iEvent);
    WCSimRootTrigger* wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

    //Searching the main track
    std::cout << "Searching the main track of event " << iEvent << " in the given file..." << "\n\n";
    WCSimRootTrack* wcsimroottrack = new WCSimRootTrack();
    int itrack;
    int ntrack = wcsimrootevent->GetNtrack_slots();
    for (itrack=0; itrack<ntrack; itrack++) {
        wcsimroottrack = dynamic_cast<WCSimRootTrack*> (wcsimrootevent->GetTracks()->At(itrack));   //Casting TObject* into WCSimRootTrack*
        if (wcsimroottrack->GetId()==1) {
            std::cout << "Main track mumber: " << itrack << "\n";
            break;
        }
    }
    if (itrack==ntrack) {
        std::cout << "Error: losing main track" << "\n";
        return;
    }

    //printing basic information
    std::cout << "Track ID: " << wcsimroottrack->GetId() << "\n" <<
        "Track ipnu (PDG code): " << wcsimroottrack->GetIpnu() << "\n" <<
        "PDG code of parent particle (expected to be 0): " << wcsimroottrack->GetParenttype() << "\n" <<
        "Track initial relativistic energy [MeV]: " << wcsimroottrack->GetE() << "\n" <<
        "Track mass [MeV/c^2]: " << wcsimroottrack->GetM() << "\n" <<
        "Starting position [cm]: (" << wcsimroottrack->GetStart(0) << ", " <<
            wcsimroottrack->GetStart(1) << ", " <<
            wcsimroottrack->GetStart(2) << ")\n" <<
        "Stopping position [cm]: (" << wcsimroottrack->GetStop(0) << ", " <<
            wcsimroottrack->GetStop(1) << ", " <<
            wcsimroottrack->GetStop(2) << ")\n" <<
        "Track initial direction [unit 3-vector]: (" << wcsimroottrack->GetDir(0) << ", " <<
            wcsimroottrack->GetDir(1) << ", " <<
            wcsimroottrack->GetDir(2) << ")\n\n";

    //printing crossing information
    int nCrossing = wcsimroottrack->GetBoundaryPoints().size(); //number of crossing
    std::cout << "Number of ID/OD crossings: " << nCrossing << "\n";
    if (nCrossing<=1) {
        std::cout << "The main track did not enter the tank.\n\n";
        return;
    }
    for (int iCrossing=0; iCrossing<nCrossing; iCrossing++) {
        std::cout << "Crossing point: " << iCrossing << "\n" << 
            "Crossing position [mm]: (" << wcsimroottrack->GetBoundaryPoints().at(iCrossing).at(0) << ", " <<
                wcsimroottrack->GetBoundaryPoints().at(iCrossing).at(1) << ", " <<
                wcsimroottrack->GetBoundaryPoints().at(iCrossing).at(2) << ")\n" <<
            "KE at crossing position [MeV]: " << wcsimroottrack->GetBoundaryKEs().at(iCrossing) << "\n" <<
            "Crossing time [ns]: " << wcsimroottrack->GetBoundaryTimes().at(iCrossing) << "\n" <<
            "Crossing type: " << wcsimroottrack->GetBoundaryTypes().at(iCrossing) << "\n\n";
    }

    //if the number of crossing point is smaller than 4, the track didn't leave the tank
    if (nCrossing<4) {
        std::cout << "The track didn't leave the tank!";
        return;
    }

    //calculate the travel distance of the main track in the tank
    double dx = wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(0) - wcsimroottrack->GetBoundaryPoints().at(0).at(0);
    double dy = wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(1) - wcsimroottrack->GetBoundaryPoints().at(0).at(1);
    double dz = wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(2) - wcsimroottrack->GetBoundaryPoints().at(0).at(2);
    std::cout << "Travel distance (estimated) of the main track in the tank [mm]: " << sqrt(dx*dx+dy*dy+dz*dz) << "\n\n";


    //Getting the geometry
    TTree* geotree = (TTree*)file->Get("wcsimGeoT");
    WCSimRootGeom* geo = new WCSimRootGeom();
    geotree->SetBranchAddress("wcsimrootgeom", &geo);
    geotree->GetEntry(0);

    //Getting the PMT hits
    int nPMTs_type0 = geo->GetWCNumPMT();
    std::vector<double> pmt_hit(nPMTs_type0, 0);    //vector storing the number of hits each pmt tubes
    int nhits = wcsimrootevent->GetNcherenkovdigihits();
    for (int i=0; i<nhits; i++) {
        WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit = (WCSimRootCherenkovDigiHit*) (wcsimrootevent->GetCherenkovDigiHits())->At(i);
        int tubeNumber = wcsimrootcherenkovdigihit->GetTubeId()-1;
        pmt_hit[tubeNumber] += wcsimrootcherenkovdigihit->GetQ();
    }

    //projection of the PMTs charge with the corresponding position on a 2D histogram
    TH2D* hist_event_display = new TH2D("Charges","Charges",250,-TMath::Pi()*max_r,TMath::Pi()*max_r,250,-max_z-2*max_r,max_z+2*max_r);
    double barrelCut = max_z-20;  //points with the absolute value of z coordinate larger than barrelCut will be treated as points on either the top or bottom of the tank 
    for (int i=0; i<nPMTs_type0; i++) {
        //rotation for event display
        double x = -(geo->GetPMT(i)).GetPosition(0);
        double y = (geo->GetPMT(i)).GetPosition(2);
        double z = (geo->GetPMT(i)).GetPosition(1);

        //barrel
        if (fabs(z)<barrelCut) {hist_event_display->Fill(-max_r*atan2(y, x), z, pmt_hit[i]);}
        //top
        else if (z>barrelCut) {hist_event_display->Fill(-y, max_z+max_r-x, pmt_hit[i]);}
        //bottom
        else {hist_event_display->Fill(-y, -max_z-max_r+x, pmt_hit[i]);}
    }
    TCanvas* c1 = new TCanvas();
    gStyle->SetOptStat(0);
    hist_event_display->Draw("colz");

    double evtx, evty;  //use as temporary variables
    //rotation for event display
    double entrance_x = -(wcsimroottrack->GetBoundaryPoints().at(0).at(0))/10;
    double entrance_y = (wcsimroottrack->GetBoundaryPoints().at(0).at(2))/10;
    double entrance_z = (wcsimroottrack->GetBoundaryPoints().at(0).at(1))/10;
    //Draw the entrance point
    //barrel
    if (fabs(entrance_z)<barrelCut) {
        evtx = -max_r*atan2(entrance_y, entrance_x);
        evty = entrance_z;
    }
    //top
    else if (entrance_z>barrelCut) {
        evtx = -entrance_y;
        evty = max_z+max_r-entrance_x;
    }
    //bottom
    else {
        evtx = -entrance_y;
        evty = -max_z-max_r+entrance_x;
    }
    TMarker m1(evtx,evty,29);
    m1.SetMarkerColor(kRed);
    m1.Draw();

    //rotation for event display
    double exit_x = -(wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(0))/10;
    double exit_y = (wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(2))/10;
    double exit_z = (wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(1))/10;
    //Draw the exit point
    //barrel
    if (fabs(exit_z)<barrelCut) {
        evtx = -max_r*atan2(exit_y, exit_x);
        evty = exit_z;
    }
    //top
    else if (exit_z>barrelCut) {
        evtx = -exit_y;
           evty = max_z+max_r-exit_x;
    }
    //bottom
    else {
        evtx = -exit_y;
        evty = -max_z-max_r+exit_x;
    }
    TMarker m2(evtx,evty,29);
    m2.SetMarkerColor(kBlack);
    m2.Draw();
    c1->SaveAs(Form("Charge_Display.pdf"));

    file->Close();
}