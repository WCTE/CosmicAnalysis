//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542
R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

//Selection option: 1 for entrance radius cutoff, 2 for exit radius cutoff
static bool selector(const int SelectOpt, double fitQunEntrance[3], const double fitQunExit[3], const double cut) {

    if (fitQunEntrance[2]!=max_z) {return true;}    //selecting the fitQun entrance point on the top of the tank
    if (fitQunExit[2]!=-max_z) {return true;}       //selecting the fitQun exit point at the bottom of the tank

    //entrance radius cutoff
    if (SelectOpt==1) {
        if ((fitQunEntrance[0]*fitQunEntrance[0]+fitQunEntrance[1]*fitQunEntrance[1])>=cut*cut) {return true;}
    }
    //exit radius cutoff
    if (SelectOpt==2) {
        if ((fitQunExit[0]*fitQunExit[0]+fitQunExit[1]*fitQunExit[1])>=cut*cut) {return true;}
    }

    return false;
}

void m_selection(const char* pname="./WCSim_fitQun_preprocess.root") {
    //starting cut positions must be smaller than the ending cut positions
    //total length (i.e. end-start) must be divisible by the interval
    double start = 0, end = 162, interval = 1;
    int bin = (end-start)/interval;     //bin of the histogram
    TH1D* entrance_diff_mean = new TH1D("entrance_diff_mean","Mean of entrance difference;cutoff;[cm]",bin,start,end);
    TH1D* exit_diff_mean = new TH1D("exit_diff_mean","Mean of exit difference;cutoff;[cm]",bin,start,end);
    TH1D* dist_diff_mean = new TH1D("dist_diff_mean","Mean of distance difference;cutoff;[cm]",bin,start,end);
    TH1D* percentage = new TH1D("percantage","Percentage of remaining;cutoff;percentage",bin,start,end);

    TFile* file = new TFile(pname, "read");
    if (!file->IsOpen()) {std::cout << "Error, could not open the preprocess file: " << pname << "\n";}
    TTree* tree = (TTree*)file->Get("WCSim_and_fitQun");
    double cut, fitQunDist, WCSimDist;
    double fitQunEntrance[3], fitQunExit[3], fitQunPos[3], WCSimEntrance[3], WCSimExit[3];
    tree->SetBranchAddress("WCSimDist", &WCSimDist);
    tree->SetBranchAddress("WCSimEntrance", WCSimEntrance);
    tree->SetBranchAddress("WCSimExit", WCSimExit);
    tree->SetBranchAddress("fitQunEntrance", fitQunEntrance);
    tree->SetBranchAddress("fitQunExit", fitQunExit);
    tree->SetBranchAddress("fitQunPos", fitQunPos);
    tree->SetBranchAddress("fitQunDist", &fitQunDist);

    int nEntries = tree->GetEntries();
    //set up histograms to calculate the mean values
    TH1D* entrance_diff = new TH1D("Entrance_diff","Entrance difference;[cm];Number of event",240,0,120);
    TH1D* exit_diff = new TH1D("Exit_diff","Exit difference;[cm];Number of event",240,0,120);
    TH1D* dist_diff = new TH1D("dist_diff","Distance difference;distance difference;Number of event",240,0,120);
    //read the file a multiple times with different cutoff
    for (cut=start;cut<end;cut+=interval) {
        int ignored_track = 0;
        for (int i=0; i<nEntries; i++) {
            tree->GetEntry(i);
            //selection of event
            if ((fitQunPos[0]*fitQunPos[0]+fitQunPos[1]*fitQunPos[1])>(max_r*max_r)) {
                ignored_track++;
                continue;
            }
            if (selector(fitQunEntrance, fitQunExit, cut)) {
                ignored_track++;
                continue;
            }

            double dx = WCSimEntrance[0] - fitQunEntrance[0];
            double dy = WCSimEntrance[1] - fitQunEntrance[1];
            double dz = WCSimEntrance[2] - fitQunEntrance[2];
            entrance_diff->Fill(sqrt(dx*dx+dy*dy+dz*dz));
            dx = WCSimExit[0] - fitQunExit[0];
            dy = WCSimExit[1] - fitQunExit[1];
            dz = WCSimExit[2] - fitQunExit[2];
            exit_diff->Fill(sqrt(dx*dx+dy*dy+dz*dz));
            dist_diff->Fill(fitQunDist-WCSimDist);
        }
        entrance_diff_mean->Fill(cut,entrance_diff->GetMean());
        entrance_diff->Reset();     //Reset histogram for the next loop
        exit_diff_mean->Fill(cut,exit_diff->GetMean());
        exit_diff->Reset();
        dist_diff_mean->Fill(cut,dist_diff->GetMean());
        dist_diff->Reset();
        percentage->Fill(cut,((double)(nEntries-ignored_track))/1000);
    }
    delete entrance_diff;
    delete exit_diff;
    delete dist_diff;

    //plotting
    TCanvas* c1 = new TCanvas();
    entrance_diff_mean->Draw("hist");
    c1->SaveAs(Form("entrance_diff_mean.pdf"));
    delete entrance_diff_mean;
    exit_diff_mean->Draw("hist");
    c1->SaveAs(Form("exit_diff_mean.pdf"));
    delete exit_diff_mean;
    dist_diff_mean->Draw("hist");
    c1->SaveAs(Form("dist_diff_mean.pdf"));
    delete dist_diff_mean;
    percentage->Draw("hist");
    c1->SaveAs(Form("percentage.pdf"));
    delete percentage;
    delete c1;
    tree->Reset();
    file->Close();
}