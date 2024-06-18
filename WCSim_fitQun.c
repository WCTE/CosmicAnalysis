//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

//function that perform extraploation
bool extrapolation(float fq1rpos[2][7][3], float fq1rdir[2][7][3], double entrance[3], double exit[3], double R=max_r, double Z=max_z) {
    double x = -fq1rpos[0][2][0], y = fq1rpos[0][2][2], z = fq1rpos[0][2][1];
    double r = sqrt(x*x+y*y);
    double x_dir = -fq1rdir[0][2][0], y_dir = fq1rdir[0][2][2], z_dir = fq1rdir[0][2][1];
    double r_dir = sqrt(x_dir*x_dir+y_dir*y_dir);
    double vt;

    //check whether the vertex position is located in the tank
    if ((fabs(z)<Z)&&(r<R)) {
        //search for exit point
        if (z_dir>0) {vt = (Z-z)/z_dir;}
        else {vt = (-Z-z)/z_dir;}
        exit[0] = x+x_dir*vt;
        exit[1] = y+y_dir*vt;
        exit[2] = z+z_dir*vt;
        //check if the particle passes through the barrel
        if ((exit[0]*exit[0]+exit[1]*exit[1])>(R*R)) {
            double cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
            vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
            exit[0] = x+x_dir*vt;
            exit[1] = y+y_dir*vt;
            exit[2] = z+z_dir*vt;
        }
        //search for entrance point
        x_dir = -x_dir;
        y_dir = -y_dir;
        z_dir = -z_dir;
        if (z_dir>0) {vt = (Z-z)/z_dir;}
        else {vt = (-Z-z)/z_dir;}
        entrance[0] = x+x_dir*vt;
        entrance[1] = y+y_dir*vt;
        entrance[2] = z+z_dir*vt;
        //check if the particle passes through the barrel
        if ((entrance[0]*entrance[0]+entrance[1]*entrance[1])>(R*R)) {
            double cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
            vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
            entrance[0] = x+x_dir*vt;
            entrance[1] = y+y_dir*vt;
            entrance[2] = z+z_dir*vt;
        }
        return false;
    }


    if (r>R) {
        double cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
        if (acos(cos_theta)>asin(R/r)) {return true;}
        if (z>Z) {
            if (z_dir>0) {return true;}
            vt = (z-Z)/z_dir;
        }
        if (z<(-Z)) {
            if (z_dir<0) {return true;}
            vt = (-Z-z)/z_dir;
        }
        if (z_dir>0) {vt = (Z-z)/z_dir;}
        else {vt = (-Z-z)/z_dir;}
        entrance[0] = x+x_dir*vt;
        entrance[1] = y+y_dir*vt;
        entrance[2] = z+z_dir*vt;
        if ((entrance[0]*entrance[0]+entrance[1]*entrance[1])>(R*R)) {
            vt = (r/r_dir)*(cos_theta-sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
            entrance[0] = x+x_dir*vt;
            entrance[1] = y+y_dir*vt;
            entrance[2] = z+z_dir*vt;
            if (fabs(entrance[2])>Z) {return true;}
        }
        x = entrance[0];
        y = entrance[1];
        z = entrance[2];
        if (z_dir>0) {vt = (Z-z)/z_dir;}
        else {vt = (-Z-z)/z_dir;}
        exit[0] = x+x_dir*vt;
        exit[1] = y+y_dir*vt;
        exit[2] = z+z_dir*vt;
        //check if the particle passes through the barrel
        if ((exit[0]*exit[0]+exit[1]*exit[1])>(R*R)) {
            cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
            vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
            exit[0] = x+x_dir*vt;
            exit[1] = y+y_dir*vt;
            exit[2] = z+z_dir*vt;
        }
        return false;
    }
    
    if (fabs(z)>Z) {
        if (z>Z) {
            if (z_dir>0) {return true;}
            vt = (z-Z)/z_dir;
        }
        if (z<(-Z)) {
            if (z_dir<0) {return true;}
            vt = (-Z-z)/z_dir;
        }
        entrance[0] = x+x_dir*vt;
        entrance[1] = y+y_dir*vt;
        entrance[2] = z+z_dir*vt;
        if ((entrance[0]*entrance[0]+entrance[1]*entrance[1])>(R*R)) {return true;}
        x = entrance[0];
        y = entrance[1];
        z = entrance[2];
        if (z_dir>0) {vt = (Z-z)/z_dir;}
        else {vt = (-Z-z)/z_dir;}
        exit[0] = x+x_dir*vt;
        exit[1] = y+y_dir*vt;
        exit[2] = z+z_dir*vt;
        //check if the particle passes through the barrel
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

void WCSim_fitQun(const char* fname="/work/kmtsui/wcte/cosmic/hk_flux/fq_mu-_*.root", const char* wname="/work/kmtsui/wcte/cosmic/hk_flux/wcsim_mu-_*.root") {

    //serve for case checking
    int f_error_track = 0, w_error_track = 0;

    //Get the files
    TChain* ft = new TChain("fiTQun");
    ft->Add(fname);
    TChain* wt = new TChain("wcsimT");
    wt->Add(wname);

    TH1D* exit_error_x = new TH1D("Difference in x", "Difference in x coordinate of the exit position;[cm];Number of event", 400, -100, 100);
    TH1D* exit_error_y = new TH1D("Difference in y", "Difference in coordinate of the exit position;[cm];Number of event", 400, -100, 100);
    TH1D* exit_error_z = new TH1D("Difference in z", "Difference in coordinate of the exit position;[cm];Number of event", 400, -100, 100);
    TH1D* entrance_error_x = new TH1D("Difference in x", "Difference in x coordinate of the entrance position;[cm];Number of event", 400, -100, 100);
    TH1D* entrance_error_y = new TH1D("Difference in y", "Difference in y coordinate of the entrance position;[cm];Number of event", 400, -100, 100);
    TH1D* entrance_error_z = new TH1D("Difference in z", "Difference in z coordinate of the entrance position;[cm];Number of event", 400, -100, 100);
    TH1D* angular_error = new TH1D("Angular difference", "Difference in direction;Radian;Number of event", 100, -0.1, 3.2);
    TH1D* exit_dist_error = new TH1D("Exit position difference", "Difference in exit position;[cm];Number of event", 400, 0, 200);
    TH1D* entrance_dist_error = new TH1D("Entrance position difference", "Difference in entrance position;[cm];Number of event", 400, 0, 200);
    TH1D* travel_dist_error = new TH1D("Travel distance difference", "Difference in travelled distance;[cm];Number of event", 1200, -300, 300);

    int nEntries = ft->GetEntries();
    std::cout << "Number of entries: " << nEntries << "\n";

    float fq1rpos[2][7][3];
    float fq1rdir[2][7][3];
    ft->SetBranchAddress("fq1rpos",fq1rpos);
    ft->SetBranchAddress("fq1rdir",fq1rdir);
    
    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    wt->SetBranchAddress("wcsimrootevent", &wcsimrootsuperevent);
    WCSimRootTrack* wcsimroottrack = new WCSimRootTrack();

    for (int i=0; i<nEntries; i++) {

        //fitQun part
        ft->GetEntry(i);
        double f_entrance[3];
        double f_exit[3];
        if (extrapolation(fq1rpos, fq1rdir, f_entrance, f_exit)) {
            f_error_track++;
            continue;
        }

        //WCSim part
        wt->GetEntry(i);
        WCSimRootTrigger* wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
        int ntrack = wcsimrootevent->GetNtrack_slots();
        int itrack;
        for (int itrack=0; itrack<ntrack; itrack++) {
            wcsimroottrack = dynamic_cast<WCSimRootTrack*> (wcsimrootevent->GetTracks()->At(itrack));   //Casting TObject* into WCSimRootTrack*
            if (wcsimroottrack->GetId()==1) {break;}
        }
        int nCrossing = wcsimroottrack->GetBoundaryPoints().size(); //number of crossing
        //ignoring tracks that did not enter or terminated in the tank
        if (nCrossing<4) {
            w_error_track++;
            continue;
        }

        double x_dir = (-wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(0) + wcsimroottrack->GetBoundaryPoints().at(0).at(0))/10;
        double y_dir = (wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(2) - wcsimroottrack->GetBoundaryPoints().at(0).at(2))/10;
        double z_dir = (wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(1) - wcsimroottrack->GetBoundaryPoints().at(0).at(1))/10;
        double w_travel_dist = sqrt(x_dir*x_dir + y_dir*y_dir + z_dir*z_dir);
        double norm = sqrt(x_dir*x_dir + y_dir*y_dir + z_dir*z_dir);
        x_dir = x_dir/norm;
        y_dir = y_dir/norm;
        z_dir = z_dir/norm;
        angular_error->Fill(acos(-x_dir*fq1rdir[0][2][0] + y_dir*fq1rdir[0][2][2] + z_dir*fq1rdir[0][2][1]));

        double dx = f_exit[0] + wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(0)/10;
        double dy = f_exit[1] - wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(2)/10;
        double dz = f_exit[2] - wcsimroottrack->GetBoundaryPoints().at(nCrossing-1).at(1)/10;
        exit_error_x->Fill(dx);
        exit_error_y->Fill(dy);
        exit_error_z->Fill(dz);
        exit_dist_error->Fill(sqrt(dx*dx+dy*dy+dz*dz));
        dx = f_entrance[0] + wcsimroottrack->GetBoundaryPoints().at(0).at(0)/10;
        dy = f_entrance[1] - wcsimroottrack->GetBoundaryPoints().at(0).at(2)/10;
        dz = f_entrance[2] - wcsimroottrack->GetBoundaryPoints().at(0).at(1)/10;
        entrance_error_x->Fill(dx);
        entrance_error_y->Fill(dy);
        entrance_error_z->Fill(dz);
        entrance_dist_error->Fill(sqrt(dx*dx+dy*dy+dz*dz));

        double f_travel_dist = sqrt((f_exit[0]-f_entrance[0])*(f_exit[0]-f_entrance[0]) + (f_exit[1]-f_entrance[1])*(f_exit[1]-f_entrance[1]) + (f_exit[2]-f_entrance[2])*(f_exit[2]-f_entrance[2]));
        travel_dist_error->Fill(f_travel_dist-w_travel_dist);
    }

    //Case checking
    std::cout << "Number of main track which did not enter the tank in WCSim files: " << w_error_track << "\n";
    std::cout << "Number of reconstructed track which did not enter the tank in fitQun files: " << f_error_track << "\n";

    TCanvas* c1 = new TCanvas();
    exit_error_x->Draw();
    c1->SaveAs(Form("Exit_diff_x.pdf"));
    delete exit_error_x;
    exit_error_y->Draw();
    c1->SaveAs(Form("Exit_diff_y.pdf"));
    delete exit_error_y;
    exit_error_z->Draw();
    c1->SaveAs(Form("Exit_diff_z.pdf"));
    delete exit_error_z;
    entrance_error_x->Draw();
    c1->SaveAs(Form("Entrance_diff_x.pdf"));
    delete entrance_error_x;
    entrance_error_y->Draw();
    c1->SaveAs(Form("Entrance_diff_y.pdf"));
    delete entrance_error_y;
    entrance_error_z->Draw();
    c1->SaveAs(Form("Entrance_diff_z.pdf"));
    delete entrance_error_z;
    angular_error->Draw();
    c1->SaveAs(Form("Angular_diff.pdf"));
    delete angular_error;
    exit_dist_error->Draw();
    c1->SaveAs(Form("Exit_diff.pdf"));
    delete exit_dist_error;
    entrance_dist_error->Draw();
    c1->SaveAs(Form("Entrance_diff.pdf"));
    delete entrance_dist_error;
    travel_dist_error->Draw();
    c1->SaveAs(Form("Travel_dist_diff.pdf"));
    delete travel_dist_error;

    ft->Reset();
    wt->Reset();
}