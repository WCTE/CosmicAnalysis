//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

void multiple_events_reader(const char* filename="/work/kmtsui/wcte/cosmic/hk_flux/wcsim_mu-_*.root") {

    //losing_track records the number of event which failed to find the main track, track_ignored records the number of main track which fails to pass through the tank
    int losing_track = 0, track_ignored = 0;

    //Get the files
    TChain* t = new TChain("wcsimT");
    t->Add(filename);

    //projection of the extrance and exit points on.2D histograms
    TH2D* entrance_display = new TH2D("Entrance_position","Entrance position;[cm];[cm]",300,-TMath::Pi()*max_r,TMath::Pi()*max_r,300,-max_z-2*max_r,max_z+2*max_r);
    TH2D* exit_display = new TH2D("Exit_position","Exit position;[cm];[cm]",300,-TMath::Pi()*max_r,TMath::Pi()*max_r,300,-max_z-2*max_r,max_z+2*max_r);
    //2D histogram of total charge against energy loss
    TH2D* Q_vs_E = new TH2D("Q_vs_E", "Total charge against energy loss;Energy dissipation [Mev];Total charge;", 3000, 0, 15000, 10000, 0, 40000);
    double barrelCut = max_z-0.1;  //points with the absolute value of z coordinate larger than barrelCut will be treated as points on either the top or bottom of the tank
    
    int nEvent = t->GetEntries();
    std::cout << "Total number of event in all the given files: " << nEvent << "\n\n";
    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    t->SetBranchAddress("wcsimrootevent", &wcsimrootsuperevent);
    WCSimRootTrack* wcsimroottrack = new WCSimRootTrack();

    for (int i=0; i<nEvent; i++) {
        //finding the main track of this event
        t->GetEntry(i);
        WCSimRootTrigger* wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
        int ntrack = wcsimrootevent->GetNtrack_slots();
        int itrack;
        for (int itrack=0; itrack<ntrack; itrack++) {
            wcsimroottrack = dynamic_cast<WCSimRootTrack*> (wcsimrootevent->GetTracks()->At(itrack));   //Casting TObject* into WCSimRootTrack*
            if (wcsimroottrack->GetId()==1) {break;}
        }
        //Search for next event if main track is missed
        if (itrack==ntrack) {
            losing_track++;
            continue;
        }
        
        int nCrossing = wcsimroottrack->GetBoundaryPoints().size(); //number of crossing
        //ignoring tracks that did not enter or terminated in the tank
        if (nCrossing<4) {
            track_ignored++;
            continue;
        }

        //Draw entrance point
        double x = -wcsimroottrack->GetBoundaryPoints().at(0).at(0)/10;
        double y = wcsimroottrack->GetBoundaryPoints().at(0).at(2)/10;
        double z = wcsimroottrack->GetBoundaryPoints().at(0).at(1)/10;
        //barrel
        if (fabs(z)<barrelCut) entrance_display->Fill(-max_r*atan2(y, x), z);
        //top
        else if (z>barrelCut) entrance_display->Fill(-y, max_z+max_r-x);
        //bottom
        else entrance_display->Fill(-y, -max_z-max_r+x);

        //Draw exit point
        x = -wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(0)/10;
        y = wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(2)/10;
        z = wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(1)/10;
        //barrel
        if (fabs(z)<barrelCut) exit_display->Fill(-max_r*atan2(y, x), z);
        //top
        else if (z>barrelCut) exit_display->Fill(-y, max_z+max_r-x);
        //bottom
        else exit_display->Fill(-y, -max_z-max_r+x);

        //Ploting total charge vs energy loss diagram
        double E_loss = wcsimroottrack->GetBoundaryKEs().at(0) - wcsimroottrack->GetBoundaryKEs().at(nCrossing-1);
        double total_charge = 0;
        int nhits = wcsimrootevent->GetNcherenkovdigihits();
        for (int i=0; i<nhits; i++) {
            WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit = (WCSimRootCherenkovDigiHit*) (wcsimrootevent->GetCherenkovDigiHits())->At(i);
            total_charge += wcsimrootcherenkovdigihit->GetQ();
        }
        Q_vs_E->Fill(E_loss, total_charge);
    }

    //Case checking
    if (losing_track>0) std::cout << "Main track not find in " << losing_track << " event\n";
    if (track_ignored>0) std::cout << "Number of track did not enter, or terminated in the tank: " << track_ignored << "\n";

    //Drawing histogram
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    entrance_display->Draw("colz");
    c1->SaveAs(Form("Entrance_Position_WCSim.pdf"));
    delete entrance_display;
    exit_display->Draw("colz");
    c1->SaveAs(Form("Exit_Position_WCSim.pdf"));
    delete exit_display;
    Q_vs_E->Draw("colz");
    c1->SaveAs(Form("Charge_vs_energy_WCSim.pdf"));

    t->Reset();
}