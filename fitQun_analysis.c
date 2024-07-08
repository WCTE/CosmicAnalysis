//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

//function for extrapolation of the reconstructed vertex, return true when the track formed by extrapolation of vertex does not across the tank
bool extrapolation(float fq1rpos[2][7][3], float fq1rdir[2][7][3], double entrance[3], double exit[3], double R=max_r, double Z=max_z) {
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

void fitQun_analysis(const char* filename = "/work/kmtsui/wcte/cosmic/hk_flux/fq_mu-_*.root") {

    //Get the files
    TChain* t = new TChain("fiTQun");
    t->Add(filename);

    int nEntries = t->GetEntries();
    std::cout << "Number of entries: " << nEntries << "\n";

    //Set up branches for fitQun data
    float fq1rpos[2][7][3];
    float fq1rdir[2][7][3];
    t->SetBranchAddress("fq1rpos",fq1rpos);
    t->SetBranchAddress("fq1rdir",fq1rdir);
    
    //projection of the extrance and exit points on.2D histograms
    TH2D* entrance_display = new TH2D("Entrance_position","Entrance position;[cm];[cm]",300,-TMath::Pi()*max_r,TMath::Pi()*max_r,300,-max_z-2*max_r,max_z+2*max_r);
    TH2D* exit_display = new TH2D("Exit_position","Exit position;[cm];[cm]",300,-TMath::Pi()*max_r,TMath::Pi()*max_r,300,-max_z-2*max_r,max_z+2*max_r);
    double barrelCut = max_z-0.01;  //points with the absolute value of z coordinate larger than barrelCut will be treated as points on either the top or bottom of the tank 

    int error_track = 0;    //error_track counts the number of tracks which do not across the tank

    for (int i=0; i<nEntries; i++) {
        t->GetEntry(i);
        double entrance[3];
        double exit[3];
        if (extrapolation(fq1rpos, fq1rdir, entrance, exit)) {
            error_track++;
            continue;
        }

        //Draw entrance point
        double x = entrance[0];
        double y = entrance[1];
        double z = entrance[2];
        //barrel
        if (fabs(z)<barrelCut) entrance_display->Fill(-max_r*atan2(y, x), z);
        //top
        else if (z>barrelCut) entrance_display->Fill(-y, max_z+max_r-x);
        //bottom
        else entrance_display->Fill(-y, -max_z-max_r+x);

        //Draw exit point
        x = exit[0];
        y = exit[1];
        z = exit[2];
        //barrel
        if (fabs(z)<barrelCut) exit_display->Fill(-max_r*atan2(y, x), z);
        //top
        else if (z>barrelCut) exit_display->Fill(-y, max_z+max_r-x);
        //bottom
        else exit_display->Fill(-y, -max_z-max_r+x);
    }

    std::cout << "Outside the tank: " << error_track << "\n";

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    entrance_display->Draw("colz");
    c1->SaveAs(Form("Entrance_Position_fitQun.pdf"));
    delete entrance_display;
    exit_display->Draw("colz");
    c1->SaveAs(Form("Exit_Position_fitQun.pdf"));
    delete exit_display;

    t->Reset();
    delete c1;
}