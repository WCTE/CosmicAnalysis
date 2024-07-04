//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

bool selector(double fitQunPos[3], double fitQunEntrance[3], double fitQunExit[3]) {
    if ((fitQunPos[0]*fitQunPos[0]+fitQunPos[1]*fitQunPos[1])<(max_r*max_r)) {return true;}
    //if (fitQunEntrance[2]!=max_z) {return true;}
    //if (fitQunExit[2]!=-max_z) {return true;}

    //entrance/exit cut
    //if ((fitQunEntrance[0]*fitQunEntrance[0]+fitQunEntrance[1]*fitQunEntrance[1])>=(120*120)) {return true;}
    //if ((fitQunExit[0]*fitQunExit[0]+fitQunExit[1]*fitQunExit[1])>=(145*145)) {return true;}

    //reconstructed distance cut
    //double dx = exit[0] - entrance[0];
    //double dy = exit[1] - entrance[1];
    //double dz = exit[2] - entrance[2];
    //if ((dx*dx+dy*dy+dz*dz)<(340*340)) {return true;}

    return false;
}

void s_selection(const char* pname="./WCSim_fitQun_preprocess.root") {
    //Get the preprocessed file
    TFile* file = new TFile(pname, "read");
    if (!file->IsOpen()){
        std::cout << "Error, could not open the preprocess file: " << pname << "\n";
        return;
    }
    TTree* tree = (TTree*)file->Get("WCSim_and_fitQun");

    double WCSimEntrance[3], WCSimExit[3], fitQunPos[3], fitQunEntrance[3], fitQunExit[3];
    double WCSimTotalQ, WCSimDist, fitQunDist, fitQunLikelihood, eta;
    tree->SetBranchAddress("fitQunPos", fitQunPos);
    tree->SetBranchAddress("WCSimEntrance", WCSimEntrance);
    tree->SetBranchAddress("WCSimExit", WCSimExit);
    tree->SetBranchAddress("fitQunEntrance", fitQunEntrance);
    tree->SetBranchAddress("fitQunExit", fitQunExit);
    tree->SetBranchAddress("WCSimDist", &WCSimDist);
    tree->SetBranchAddress("fitQunDist", &fitQunDist);
    //tree->SetBranchAddress("fitQunLikelihood", &fitQunLikelihood);
    //tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("WCSimTotalQ", &WCSimTotalQ);

    TH1D* entrance_diff = new TH1D("Entrance_diff","Entrance difference;[cm];Number of event",240,0,120);
    TH1D* exit_diff = new TH1D("Exit_diff","Exit difference;[cm];Number of event",240,0,120);
    TH1D* dist_diff = new TH1D("dist_diff","Distance difference;distance difference [cm];Number of event",480,-120,120);
    TH2D* QvsRL = new TH2D("QvsRL","Total charge against reconstructed distance;reconstructed distance [cm];total charge",45,0,435,110,0,15000);
    TH2D* QvsTL = new TH2D("QvsTL","Total charge against true distance;true distance [cm];total charge",80,0,435,100,0,15000);
    TH1D* Q_over_trueL = new TH1D("Q_over_trueL","Total charge over true distance;charge over true distance;Number of events",100,0,50);
    TH1D* Q_over_reconstructedL = new TH1D("Q_over_reconstructedL","Total charge over reconstructed distance;charge over reconstructed distance;Number of events",100,0,50);
    
    int nEntries = tree->GetEntries();
    std::cout << "Number of entries: " << nEntries << "\n";

    //serve for case checking
    int ignored_track = 0;
    for (int i=0; i<nEntries; i++) {
        tree->GetEntry(i);
        //selection of event
        if (selector(fitQunPos, fitQunEntrance, fitQunExit)) {
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
        QvsRL->Fill(fitQunDist,WCSimTotalQ);
        QvsTL->Fill(WCSimDist,WCSimTotalQ);
        Q_over_trueL->Fill(WCSimTotalQ/WCSimDist);
        Q_over_reconstructedL->Fill(WCSimTotalQ/fitQunDist);
    }

    //Case checking
    std::cout << "Number of reconstructed tracks ignored: " << ignored_track << "\n";

    TCanvas* c1 = new TCanvas();
    entrance_diff->Draw();
    c1->SaveAs(Form("entrance_diff.pdf"));
    delete entrance_diff;
    exit_diff->Draw();
    c1->SaveAs(Form("exit_diff.pdf"));
    delete exit_diff;
    dist_diff->Draw();
    c1->SaveAs(Form("dist_diff.pdf"));
    delete dist_diff;
    Q_over_trueL->Draw();
    c1->SaveAs(Form("Q_over_trueL.pdf"));
    Q_over_reconstructedL->Draw();
    c1->SaveAs(Form("Q_over_reconstructedL.pdf"));

    gStyle->SetOptStat(0);
    Q_over_trueL->SetStats(0);
    Q_over_trueL->SetLineColor(kRed);
    Q_over_trueL->SetTitle("Total charge over distance");
    Q_over_trueL->SetXTitle("charge over distance");
    Q_over_trueL->Draw();
    Q_over_reconstructedL->SetLineColor(kBlue);
    Q_over_reconstructedL->Draw("same");
    TLegend* legend = new TLegend(0.75,0.75,0.9,0.9);
    legend->AddEntry(Q_over_trueL, "true", "l");
    legend->AddEntry(Q_over_reconstructedL, "reconstructed", "l");
    legend->SetTextSize(0.04);
    legend->Draw();
    c1->SaveAs(Form("Q_over_dist.pdf"));
    delete Q_over_trueL;
    delete Q_over_reconstructedL;
    delete legend;

    QvsRL->Draw("Colz");
    c1->SaveAs(Form("QvsRL.pdf"));
    delete QvsRL;
    QvsTL->Draw("Colz");
    c1->SaveAs(Form("QvsTL.pdf"));
    delete QvsTL;
    delete c1;

    tree->Reset();
    file->Close();
}