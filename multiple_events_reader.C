//max_r and max_z are obtained from find_max_r_and_z.C
#define max_r 190.446
#define max_z 156.013

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

int multiple_events_reader(const char* filename="/work/kmtsui/wcte/cosmic/hk_flux/wcsim_mu-_*.root") {

    //serve for case checking
    int losing_track = 0, track_ignored = 0;

    //Get the files
    TChain* t = new TChain("wcsimT");
    t->Add(filename);

    //setup of histogram
    TH2D* entrance_display = new TH2D("Entrance_position","Entrance position;[cm];[cm]",300,-TMath::Pi()*max_r,TMath::Pi()*max_r,300,-max_z-2*max_r,max_z+2*max_r);
    TH2D* exit_display = new TH2D("Exit_position","Exit position;[cm];[cm]",300,-TMath::Pi()*max_r,TMath::Pi()*max_r,300,-max_z-2*max_r+5,max_z+2*max_r+5);
    TH1D* Energy_dissipation = new TH1D("Energy_dissipation", "Energy dissipation distribution;[MeV/cm];Number of event", 200, 0, 10);
    TH2D* Q_vs_E = new TH2D("Q_vs_E", "Total charge against energy loss;Energy dissipation [Mev];Total charge;", 3000, 0, 15000, 10000, 0, 40000);
    double barrelCut = max_z-0.1;
    
    int nEvent = t->GetEntries();
    std::cout << "Total number of event in all the given files: " << nEvent << "\n\n";
    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    t->SetBranchAddress("wcsimrootevent", &wcsimrootsuperevent);
    WCSimRootTrack* wcsimroottrack = new WCSimRootTrack();

    for (int i=0; i<nEvent; i++) {

        std::cout << "\rReading event " << i << "/" << nEvent;

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
        if (nCrossing<3) {
            track_ignored++;
            continue;
        }

        //Determine the location of the entrance and exit point
        int entrance_boundary_index, exit_boundary_index;
        if (nCrossing<5) {
            entrance_boundary_index = 1;
            exit_boundary_index = 2;
        }
        else {
            entrance_boundary_index = 0;
            exit_boundary_index = 4;
        }

        //Draw entrance point
        double x = -wcsimroottrack->GetBoundaryPoints().at(entrance_boundary_index).at(0)/10;
        double y = wcsimroottrack->GetBoundaryPoints().at(entrance_boundary_index).at(2)/10;
        double z = wcsimroottrack->GetBoundaryPoints().at(entrance_boundary_index).at(1)/10;
        //barrel
        if (fabs(z)<barrelCut) entrance_display->Fill(-max_r*atan2(y, x), z, 1);
        //top
        else if (z>barrelCut) entrance_display->Fill(-y, max_z+max_r-x, 1);
        //bottom
        else entrance_display->Fill(-y, -max_z-max_r+x, 1);

        //Draw exit point
        x = -wcsimroottrack->GetBoundaryPoints().at(exit_boundary_index).at(0)/10;
        y = wcsimroottrack->GetBoundaryPoints().at(exit_boundary_index).at(2)/10;
        z = wcsimroottrack->GetBoundaryPoints().at(exit_boundary_index).at(1)/10;
        //barrel
        if (fabs(z)<barrelCut) exit_display->Fill(-max_r*atan2(y, x), z, 1);
        //top
        else if (z>barrelCut) exit_display->Fill(-y, max_z+max_r-x, 1);
        //bottom
        else exit_display->Fill(-y, -max_z-max_r+x, 1);

        //Ploting total charge vs energy loss diagram
        double E_loss = wcsimroottrack->GetBoundaryKEs().at(entrance_boundary_index) - wcsimroottrack->GetBoundaryKEs().at(exit_boundary_index);
        double total_charge = 0;
        int nhits = wcsimrootevent->GetNcherenkovdigihits();
        for (int i=0; i<nhits; i++) {
            WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit = (WCSimRootCherenkovDigiHit*) (wcsimrootevent->GetCherenkovDigiHits())->At(i);
            total_charge += wcsimrootcherenkovdigihit->GetQ();
        }
        Q_vs_E->Fill(E_loss, total_charge, 1);
    }

    //Case checking
    if (losing_track>0) std::cout << "Main track not find in " << losing_track << " event\n";
    if (track_ignored>0) std::cout << "Number of track did not enter, or terminated in the tank: " << track_ignored << "\n";

    //Drawing histogram
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    entrance_display->Draw("colz");
    c1->SaveAs(Form("Entrance_Position.pdf"));

    TCanvas* c2 = new TCanvas();
    exit_display->Draw("colz");
    c2->SaveAs(Form("Exit_Position.pdf"));

    TCanvas* c3 = new TCanvas();
    Q_vs_E->Draw("colz");
    c3->SaveAs(Form("Charge_vs_energy.pdf"));

    return 0;
}