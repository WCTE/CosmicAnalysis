//max_r and max_z is taken to be the outest boundary of the tank
#define max_r 162.67
#define max_z 143.542

R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

template<class T> vector<double> GetData(const char* pname) {
    //Get the preprocessed file
    TFile* file = new TFile(pname, "read");
    if (!file->IsOpen()) {std::cout << "Error, could not open the preprocess file: " << pname << "\n";}
    TTree* tree = (TTree*)file->Get("WCSim_and_fitQun");

    TH1D* entrance_diff = new TH1D("Entrance_diff","Entrance difference;[cm];Number of event",240,0,120);
    TH1D* exit_diff = new TH1D("Exit_diff","Exit difference;[cm];Number of event",240,0,120);
    TH1D* dist_diff = new TH1D("dist_diff","Distance difference;distance difference;Number of event",240,0,120);
    
    double WCSimDist, WCSimEntrance[3], WCSimExit[3];
    tree->SetBranchAddress("WCSimDist", &WCSimDist);
    tree->SetBranchAddress("WCSimEntrance", WCSimEntrance);
    tree->SetBranchAddress("WCSimExit", WCSimExit);
    tree->SetBranchAddress("fitQunEntrance", T::fitQunEntrance);
    tree->SetBranchAddress("fitQunExit", T::fitQunExit);
    tree->SetBranchAddress("fitQunDist", &T::fitQunDist);
    tree->SetBranchAddress("fitQunLikelihood", &T::fitQunLikelihood);

    int nEntries = tree->GetEntries();
    int ignored_track = 0;
    vector<double> data;
    for (int i=0; i<nEntries; i++) {
        tree->GetEntry(i);
        //selection of event
        if (T::selector()) {
            ignored_track++;
            continue;
        }

        double dx = WCSimEntrance[0] - T::fitQunEntrance[0];
        double dy = WCSimEntrance[1] - T::fitQunEntrance[1];
        double dz = WCSimEntrance[2] - T::fitQunEntrance[2];
        double entrance_delta = sqrt(dx*dx+dy*dy+dz*dz);
        entrance_diff->Fill(entrance_delta);
        dx = WCSimExit[0] - T::fitQunExit[0];
        dy = WCSimExit[1] - T::fitQunExit[1];
        dz = WCSimExit[2] - T::fitQunExit[2];
        double exit_delta = sqrt(dx*dx+dy*dy+dz*dz);
        exit_diff->Fill(exit_delta);
        dist_diff->Fill(fabs(WCSimDist-T::fitQunDist));
    }
    data.push_back(entrance_diff->GetMean());
    data.push_back(entrance_diff->GetStdDev());
    delete entrance_diff;
    data.push_back(exit_diff->GetMean());
    data.push_back(exit_diff->GetStdDev());
    delete exit_diff;
    data.push_back(dist_diff->GetMean());
    data.push_back(dist_diff->GetStdDev());
    delete dist_diff;
    data.push_back(((double)(nEntries-ignored_track))/1000);

    tree->Reset();
    file->Close();
    return data;
}
    
class selector_2D {
public:
    static double cut, cut2, fitQunDist, fitQunLikelihood;
    static double fitQunEntrance[3], fitQunExit[3];
    selector_2D(const char* pname, double start, double end, double interval, double start2, double end2, double interval2) {
        int bin = (end-start)/interval, bin2 = (end2-start2)/interval2;
        TH2D* entrance_diff_mean = new TH2D("entrance_diff_mean","Mean of entrance difference;cutoff 1;cutoff 2",bin,start,end,bin2,start2,end2);
        TH2D* entrance_diff_sd = new TH2D("entrance_diff_sd","SD of entrance difference;cutoff 1;cutoff 2",bin,start,end,bin2,start2,end2);
        TH2D* exit_diff_mean = new TH2D("exit_diff_mean","Mean of exit difference;cutoff 1;cutoff 2",bin,start,end,bin2,start2,end2);
        TH2D* exit_diff_sd = new TH2D("exit_diff_sd","SD of exit difference;cutoff 1;cutoff 2",bin,start,end,bin2,start2,end2);
        TH2D* dist_diff_mean = new TH2D("dist_diff_mean","Mean of distance difference;cutoff 1;cutoff 2",bin,start,end,bin2,start2,end2);
        TH2D* dist_diff_sd = new TH2D("dist_diff_sd","SD of distance difference;cutoff 1;cutoff 2",bin,start,end,bin2,start2,end2);
        TH2D* percentage = new TH2D("percantage","Percentage of remaining;cutoff 1;cutoff 2",bin,start,end,bin2,start2,end2);
        
        for (cut=start;cut<end;cut+=interval) {
            for (cut2=start2;cut2<end2;cut2+=interval2) {
                vector<double> data = GetData<selector_2D> (pname);
                entrance_diff_mean->Fill(cut,cut2,data[0]);
                entrance_diff_sd->Fill(cut,cut2,data[1]);
                exit_diff_mean->Fill(cut,cut2,data[2]);
                exit_diff_sd->Fill(cut,cut2,data[3]);
                dist_diff_mean->Fill(cut,cut2,data[4]);
                dist_diff_sd->Fill(cut,cut2,data[5]);
                percentage->Fill(cut,cut2,data[6]);
            }
        }

        TCanvas* c1 = new TCanvas();
        gStyle->SetOptStat(0);
        entrance_diff_mean->Draw("Colz");
        c1->SaveAs(Form("entrance_diff_mean.pdf"));
        delete entrance_diff_mean;
        entrance_diff_sd->Draw("Colz");
        c1->SaveAs(Form("entrance_diff_sd.pdf"));
        delete entrance_diff_sd;
        exit_diff_mean->Draw("Colz");
        c1->SaveAs(Form("exit_diff_mean.pdf"));
        delete exit_diff_mean;
        exit_diff_sd->Draw("Colz");
        c1->SaveAs(Form("exit_diff_sd.pdf"));
        delete exit_diff_sd;
        dist_diff_mean->Draw("Colz");
        c1->SaveAs(Form("dist_diff_mean.pdf"));
        delete dist_diff_mean;
        dist_diff_sd->Draw("Colz");
        c1->SaveAs(Form("dist_diff_sd.pdf"));
        delete dist_diff_sd;
        percentage->Draw("Colz");
        c1->SaveAs(Form("percentage.pdf"));
        delete percentage;
    }
    static bool selector() {
        if ((fitQunEntrance[0]*fitQunEntrance[0]+fitQunEntrance[1]*fitQunEntrance[1])>=cut*cut) {return true;}
        //if ((fitQunExit[0]*fitQunExit[0]+fitQunExit[1]*fitQunExit[1])>=cut2*cut2) {return true;}
        if (fitQunExit[2]>=cut2) {return true;}
        return false;
    }
};
double selector_2D::cut, selector_2D::cut2, selector_2D::fitQunDist, selector_2D::fitQunLikelihood;
double selector_2D::fitQunEntrance[3], selector_2D::fitQunExit[3];

class selector_1D {
public:
    static double cut, fitQunDist, fitQunLikelihood;
    static double fitQunEntrance[3], fitQunExit[3];
    selector_1D(const char* pname, double start, double end, double interval) {
        int bin = (end-start)/interval;
        TH1D* entrance_diff_mean = new TH1D("entrance_diff_mean","Mean of entrance difference;cutoff;[cm]",bin,start,end);
        TH1D* entrance_diff_sd = new TH1D("entrance_diff_sd","SD of entrance difference;cutoff;[cm]",bin,start,end);
        TH1D* exit_diff_mean = new TH1D("exit_diff_mean","Mean of exit difference;cutoff;[cm]",bin,start,end);
        TH1D* exit_diff_sd = new TH1D("exit_diff_sd","SD of exit difference;cutoff;[cm]",bin,start,end);
        TH1D* dist_diff_mean = new TH1D("dist_diff_mean","Mean of distance difference;cutoff;[cm]",bin,start,end);
        TH1D* dist_diff_sd = new TH1D("dist_diff_sd","SD of distance difference;cutoff;[cm]",bin,start,end);
        TH1D* percentage = new TH1D("percantage","Percentage of remaining;cutoff;percentage",bin,start,end);

        for (cut=start;cut<end;cut+=interval) {
            vector<double> data = GetData<selector_1D> (pname);
            entrance_diff_mean->Fill(cut,data[0]);
            entrance_diff_sd->Fill(cut,data[1]);
            exit_diff_mean->Fill(cut,data[2]);
            exit_diff_sd->Fill(cut,data[3]);
            dist_diff_mean->Fill(cut,data[4]);
            dist_diff_sd->Fill(cut,data[5]);
            percentage->Fill(cut,data[6]);
        }

        TCanvas* c1 = new TCanvas();
        gStyle->SetOptStat(0);
        entrance_diff_mean->Draw("hist");
        c1->SaveAs(Form("entrance_diff_mean.pdf"));
        delete entrance_diff_mean;
        entrance_diff_sd->Draw("hist");
        c1->SaveAs(Form("entrance_diff_sd.pdf"));
        delete entrance_diff_sd;
        exit_diff_mean->Draw("hist");
        c1->SaveAs(Form("exit_diff_mean.pdf"));
        delete exit_diff_mean;
        exit_diff_sd->Draw("hist");
        c1->SaveAs(Form("exit_diff_sd.pdf"));
        delete exit_diff_sd;
        dist_diff_mean->Draw("hist");
        c1->SaveAs(Form("dist_diff_mean.pdf"));
        delete dist_diff_mean;
        dist_diff_sd->Draw("hist");
        c1->SaveAs(Form("dist_diff_sd.pdf"));
        delete dist_diff_sd;
        percentage->Draw("hist");
        c1->SaveAs(Form("percentage.pdf"));
        delete percentage;
    }

    static bool selector() {
        //if (fitQunLikelihood>=cut) {return true;}
        if (fitQunEntrance[2]==max_z) {return true;}
        if (fitQunExit[2]==-max_z) {return true;}
        if (fitQunDist<=cut) {return true;}
        //if (sqrt(fitQunEntrance[0]*fitQunEntrance[0]+fitQunEntrance[1]*fitQunEntrance[1])>=cut) {return true;}
        //if (sqrt(fitQunExit[0]*fitQunExit[0]+fitQunExit[1]*fitQunExit[1])>=cut) {return true;}
        //if (fitQunEntrance[2]>=cut) {return true;}
        //if (fitQunExit[2]>=cut) {return true;}
        return false;
    }
};
double selector_1D::cut, selector_1D::fitQunDist, selector_1D::fitQunLikelihood;
double selector_1D::fitQunEntrance[3], selector_1D::fitQunExit[3];

void m_selection(const char* pname="./WCSim_fitQun_preprocess.root") {
    //starting cut positions must be smaller than the ending cut positions
    //total length (i.e. end-start) must be divisible by the interval
    double start = 0, end = 350, interval = 1;
    //double start2 = -143, end2 = 143, interval2 = 2;
    selector_1D* ptr = new selector_1D(pname, start, end, interval);
    //selector_2D* ptr = new selector_2D(pname, start, end, interval, start2, end2, interval2);
    delete ptr;
}