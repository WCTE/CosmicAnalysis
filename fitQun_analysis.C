//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

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
        x_dir=-x_dir, y_dir=-y_dir, z_dir=-z_dir;
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

    //check whether the vertex position is on top/at the bottom of the tank, and perform forward extrapolation
    if (z_dir>0) {vt = (-Z-z)/z_dir;}
    else {vt = (Z-z)/z_dir;}
    double cos_theta;
    if ((fabs(z)>Z)&&(vt>0)) {
        entrance[0] = x+x_dir*vt;
        entrance[1] = y+y_dir*vt;
        entrance[2] = z+z_dir*vt;
        if ((entrance[0]*entrance[0]+entrance[1]*entrance[1])>(R*R)) {return true;} //particle didn't touch the tank
        x = entrance[0];
        y = entrance[1]; 
        z = entrance[2]; 
        r = sqrt(x*x+y*y);
        cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
        vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
        exit[0] = x+x_dir*vt;
        exit[1] = y+y_dir*vt;
        exit[2] = z+z_dir*vt;
        return false;
    }

    //check whether the vertex position is within the tank z coordinate boundary, and perform forward extrapolation
    cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
    if ((fabs(z)<Z)&&(acos(fabs(cos_theta))<asin(R/r))) {
        vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
        entrance[0] = x+x_dir*vt;
        entrance[1] = y+y_dir*vt;
        entrance[2] = z+z_dir*vt;
        if (fabs(entrance[2])>Z) {return true;} //particle didn't touch the tank
        x = entrance[0];
        y = entrance[1];
        z = entrance[2];
        if (z_dir>0) {vt = (Z-z)/z_dir;}
        else {vt = (-Z-z)/z_dir;}
        exit[0] = x+x_dir*vt;
        exit[1] = y+y_dir*vt;
        exit[2] = z+z_dir*vt;
        if ((exit[0]*exit[0]+exit[1]*exit[1])>(R*R)) {
            cos_theta = -(x_dir*x+y_dir*y)/(r_dir*R);
            vt = 2*R*cos_theta/r_dir;
            exit[0] = x+x_dir*vt;
            exit[1] = y+y_dir*vt;
            exit[2] = z+z_dir*vt;
        }
        return false;
    }

    x_dir=-x_dir;
    y_dir=-y_dir;
    z_dir=-z_dir;
    //check whether the vertex position is on top/at the bottom of the tank, and perform backward extrapolation
    if (z_dir>0) {vt = (-Z-z)/z_dir;}
    else {vt = (Z-z)/z_dir;}
    if ((fabs(z)>Z)&&(vt>0)) {
        exit[0] = x+x_dir*vt; 
        exit[1] = y+y_dir*vt; 
        exit[2] = z+z_dir*vt;
        if ((exit[0]*exit[0]+exit[1]*exit[1])>(R*R)) {return true;} //particle didn't touch the tank
        x = exit[0]; 
        y = exit[1]; 
        z = exit[2]; 
        r = sqrt(x*x+y*y);
        cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
        vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
        entrance[0] = x+x_dir*vt; 
        entrance[1] = y+y_dir*vt; 
        entrance[2] = z+z_dir*vt;
        return false;
    }

    //check whether the vertex position is within the tank z coordinate boundary, and perform backward extrapolation
    cos_theta = -(x_dir*x+y_dir*y)/(r_dir*r);
    if ((fabs(z)<Z)&&(acos(fabs(cos_theta))<asin(R/r))) {
        vt = (r/r_dir)*(cos_theta+sqrt(cos_theta*cos_theta+(R*R)/(r*r)-1));
        exit[0] = x+x_dir*vt; 
        exit[1] = y+y_dir*vt; 
        exit[2] = z+z_dir*vt;
        if (fabs(exit[2])>Z) {return true;} //particle didn't touch the tank
        x = exit[0];
        y = exit[1];
        z = exit[2];
        if (z_dir>0) {vt = (Z-z)/z_dir;}
        else {vt = (-Z-z)/z_dir;}
        entrance[0] = x+x_dir*vt; 
        entrance[1] = y+y_dir*vt; 
        entrance[2] = z+z_dir*vt;
        //check if the particle passes through the barrel
        if ((entrance[0]*entrance[0]+entrance[1]*entrance[1])>(R*R)) {
            cos_theta = -(x_dir*x+y_dir*y)/(r_dir*R);
            vt = 2*R*cos_theta/r_dir;
            entrance[0] = x+x_dir*vt; 
            entrance[1] = y+y_dir*vt; 
            entrance[2] = z+z_dir*vt;
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

    float fq1rpos[2][7][3];
    float fq1rdir[2][7][3];
    t->SetBranchAddress("fq1rpos",fq1rpos);
    t->SetBranchAddress("fq1rdir",fq1rdir);
    
    TH2D* entrance_display = new TH2D("Entrance_position","Entrance position;[cm];[cm]",300,-TMath::Pi()*max_r,TMath::Pi()*max_r,300,-max_z-2*max_r,max_z+2*max_r);
    TH2D* exit_display = new TH2D("Exit_position","Exit position;[cm];[cm]",300,-TMath::Pi()*max_r,TMath::Pi()*max_r,300,-max_z-2*max_r+5,max_z+2*max_r+5);
    double barrelCut = max_z-0.01;

    int error_track = 0;

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

    exit_display->Draw("colz");
    c1->SaveAs(Form("Exit_Position_fitQun.pdf"));

    t->Reset();
}