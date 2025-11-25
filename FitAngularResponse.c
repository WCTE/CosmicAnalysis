const double WaterRefractiveIndex=1.33;
const double SpeedOfLight=29.9792458; // cm/ns
const double WaterCherenkovAng=acos(1./WaterRefractiveIndex);

#include "json.hpp"
using json = nlohmann::json;

// Mapping of mPMT slot ID and type
std::vector<int> mpmt_type;

void Read_mPMT_map()
{    std::string geometry_json = "wcte_v11_20250513.json";
    // json output from Geometry package
    std::ifstream ifs(geometry_json.c_str());
    json data = json::parse(ifs);
    mpmt_type.clear();
    for (int i=0;i<106;i++) // read mPMT type from json
    {
        auto s = data.at("mpmt_type").at(i);
        if (s=="Ex-situ_W") mpmt_type.emplace_back(0);
        else if (s=="Ex-situ_T") mpmt_type.emplace_back(1);
        else if (s=="In-situ_W") mpmt_type.emplace_back(2);
        else if (s=="In-situ_T") mpmt_type.emplace_back(3);
        else mpmt_type.emplace_back(-1);
    }
}

// Calculate expected time of flight and photon angle/distance per PMT
std::vector<double> Calc_TOF_Costh_Dist(TVector3 vtx, TVector3 dir, TVector3 pmtpos, TVector3 pmtdir)
{
    TVector3 Vtx_to_PMT = pmtpos - vtx;
    double Vtx_to_PMT_Theta = Vtx_to_PMT.Angle(dir);
    
    if (Vtx_to_PMT_Theta>WaterCherenkovAng) return std::vector<double> (3,0);
    
    double muon_track_length = sin(WaterCherenkovAng - Vtx_to_PMT_Theta)*Vtx_to_PMT.Mag()/sin(TMath::Pi()-WaterCherenkovAng);
    TVector3 photon_vtx = vtx+muon_track_length*dir;
    TVector3 photon_to_PMT = pmtpos-photon_vtx;
    double photon_costh = -pmtdir.Dot(photon_to_PMT.Unit());
    double tof = (muon_track_length + photon_to_PMT.Mag()*WaterRefractiveIndex)/SpeedOfLight;

    // std::cout<<muon_track_length/SpeedOfLight<<" "<<photon_to_PMT.Mag()*WaterRefractiveIndex/SpeedOfLight<<" "<<Vtx_to_PMT.Mag()/SpeedOfLight*WaterRefractiveIndex-tof<<std::endl;

    return std::vector<double>{tof,photon_costh,photon_to_PMT.Mag()};
}

void Fit(
    const char * fname= "muon_MC_selected.root",
    const char * outname= "angular_MC.root",
    const char * suffix=  "_MC",
    const char * mask_file_name="maskFile_1766.txt"
)
{
    std::string geofile = "geofile_NuPRISMBeamTest_16cShort_mPMT.txt";
    ifstream infile(geofile.c_str());
    // skip N header lines
    std::string line;
    for (int i=0; i<5; i++) std::getline(infile, line);

    int nPMTs = 1843;
    int cab_id, mpmtid, pmtid, cycloc;
    double x, y, z, dx, dy, dz;
    std::vector<int> slot_id(nPMTs);
    std::vector<int> pos_id(nPMTs);
    std::vector<TVector3> pmt_pos(nPMTs);
    std::vector<TVector3> pmt_dir(nPMTs);
    while ( infile >> cab_id >> mpmtid >> pmtid >> x >> y >> z >> dx >> dy >> dz >> cycloc )
    {
        slot_id[cab_id-1] = mpmtid;
        pos_id[cab_id-1] = pmtid-1;
        pmt_pos[cab_id-1] = TVector3(x,y,z);
        pmt_dir[cab_id-1] = TVector3(dx,dy,dz);
    }

    ifstream maskFile(mask_file_name);
    std::vector<bool> mask_list(nPMTs,false);
    while (maskFile >> cab_id)
    {
        mask_list[cab_id] = true;
    }

    TFile* f = new TFile(fname);
    TTree* t = (TTree*)f->Get("muon");
    float vertex[4], direction[3], momentum, mu_e_llh, entrance_pos[4], exit_pos[4], mu_nll;
    float truth_entrance_pos[4], truth_exit_pos[4];
    std::vector<int> * pmt_id = 0;
    std::vector<float> * pmt_Q = 0;
    std::vector<double> * pmt_T = 0;
    t->SetBranchAddress("vertex",&vertex);
    t->SetBranchAddress("direction",&direction);
    t->SetBranchAddress("entrance_pos",&entrance_pos);
    t->SetBranchAddress("exit_pos",&exit_pos);
    t->SetBranchAddress("truth_entrance_pos",&truth_entrance_pos);
    t->SetBranchAddress("truth_exit_pos",&truth_exit_pos);
    t->SetBranchAddress("pmt_id",&pmt_id);
    t->SetBranchAddress("pmt_Q",&pmt_Q);
    t->SetBranchAddress("pmt_T",&pmt_T);
    t->SetBranchAddress("mu_e_llh",&mu_e_llh);
    t->SetBranchAddress("mu_nll",&mu_nll);

    TFile* fout = new TFile(outname,"RECREATE");
    TH1D* hist_time_residual = new TH1D("hist_time_residual","hist_time_residual",1000,-10,10);
    std::vector<TH1D*> hist_time_residual_mpmt(106);
    for (int i=0;i<106;i++)
    {
        hist_time_residual_mpmt[i] = new TH1D(Form("hist_time_residual_%i",i),Form("hist_time_residual_%i",i),1000,-50,50); // per mPMT
    }
    // Ratio of charges with (time residual < -2 ns)
    TH1D* hist_bad_hit_ratio = new TH1D("hist_bad_hit_ratio","hist_bad_hit_ratio",100,0,0.1);
    TH2D* hist_end_cap_1 = new TH2D("hist_end_cap_1","hist_end_cap_1",100,-200,200,100,-200,200);
    TH2D* hist_end_cap_2 = new TH2D("hist_end_cap_2","hist_end_cap_2",100,-200,200,100,-200,200);
    TH2D* hist_end_cap_3 = new TH2D("hist_end_cap_3","hist_end_cap_3",100,-200,200,100,-200,200);

    std::vector<TH2D*> hist_QR_Costh(4); // Histogram of PMT_Q*Photon_distance vs. Photon_Costh per mPMT type 
    std::vector<TH1D*> hist_QR_Costh_1D(4); // Mean of QR(costh) per mPMT type 
    for (int i=0;i<4;i++)
    {
        hist_QR_Costh[i] = new TH2D(Form("hist_QR_Costh_%i",i),Form("hist_QR_Costh_%i",i),100,0,1,5000,0,10000); 
        hist_QR_Costh_1D[i] = new TH1D(Form("hist_QR_Costh_1D_%i",i),Form("hist_QR_Costh_1D_%i",i),100,0,1);
    }

    std::vector<TH2D*> hist_QR_Costh_mPMT(106); // per mPMT slot ID
    std::vector<TH1D*> hist_QR_Costh_mPMT_1D(106);
    for (int i=0;i<106;i++)
    {
        hist_QR_Costh_mPMT[i] = new TH2D(Form("hist_QR_Costh_Slot_%i",i),Form("hist_QR_Costh_Slot_%i",i),20,0,1,500,0,10000); 
        hist_QR_Costh_mPMT_1D[i] = new TH1D(Form("hist_QR_Costh_1D_Slot_%i",i),Form("hist_QR_Costh_1D_Slot_%i",i),20,0,1);
    }

    int nev = t->GetEntries();
    int count1pc = nev/100;
    if (count1pc==0) count1pc=1;
    for (int i=0;i<nev;i++)
    {
        if (i%(count1pc)==0) std::cout<<"Running "<<i<<"-th event of total "<<nev<<" events"<<std::endl;
        t->GetEntry(i);

        int nhit = (*pmt_id).size();
        std::vector<double> tof(nPMTs,-1);
        std::vector<double> costh(nPMTs,-1);
        std::vector<double> dist(nPMTs,-1);
        TVector3 vtx(entrance_pos[0],entrance_pos[1],entrance_pos[2]);
        TVector3 dir(exit_pos[0]-entrance_pos[0],exit_pos[1]-entrance_pos[1],exit_pos[2]-entrance_pos[2]);
        dir = dir.Unit();

        // Calculate expected time of flight and photon angle/distance per PMT
        for (int j=0;j<nPMTs;j++)
        {
            if (mask_list[j]) continue;
            std::vector<double> vec = Calc_TOF_Costh_Dist(vtx,dir,pmt_pos[j],pmt_dir[j]);
            tof[j] = vec[0];
            costh[j] = vec[1];
            dist[j] = vec[2];
        }
        std::vector<double> QR_vec(nPMTs,0);
        
        double bad_hit=0;
        double good_hit=0;
        for (int j=0;j<nhit;j++)
        {
            cab_id = (*pmt_id)[j];
            double Q = (*pmt_Q)[j];
            double T = (*pmt_T)[j];

            if (tof[cab_id]>0)
            {
                // std::cout<<cab_id<<" "<<Q<<" "<<T<<"\n";
                double t_res = T-entrance_pos[3]-tof[cab_id];
                hist_time_residual->Fill(t_res,Q);
                mpmtid = slot_id[cab_id];
                hist_time_residual_mpmt[mpmtid]->Fill(t_res,Q);

                int mpmt_type_id = mpmt_type[mpmtid];

                // only count close-enough hits
                if (t_res>-2 && t_res<2)
                {
                    QR_vec[cab_id] += Q*dist[cab_id];
                }

                if (t_res<-2) 
                {
                    bad_hit += Q ;
                }
                else
                {
                    good_hit += Q;
                }
            }
        }

        hist_bad_hit_ratio->Fill(bad_hit/(bad_hit+good_hit));

        if (bad_hit/(bad_hit+good_hit)<0.008) 
        {
            hist_end_cap_1->Fill(exit_pos[0],exit_pos[2]);

            for (int j=0;j<nPMTs;j++)
            {
                if (tof[j]<=0) continue;
                if (dist[j]<100) continue; // only consider far enough hits in space
                mpmtid = slot_id[j];
                int mpmt_type_id = mpmt_type[mpmtid];
                hist_QR_Costh[mpmt_type_id]->Fill(costh[j],QR_vec[j]);
                hist_QR_Costh_mPMT[mpmtid]->Fill(costh[j],QR_vec[j]);
            }
        }
        else if (bad_hit/(bad_hit+good_hit)<0.02) hist_end_cap_2->Fill(exit_pos[0],exit_pos[2]);
        else hist_end_cap_3->Fill(exit_pos[0],exit_pos[2]);
    }

    hist_time_residual->Draw("hist");
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SaveAs(Form("time_residual%s.pdf",suffix));

    TH1D* hist_time_residual_fit = new TH1D("hist_time_residual_fit", "hist_time_residual_fit", 106, 0, 106);
    for (int i=0;i<106;i++)
    {
        if (hist_time_residual_mpmt[i]->GetEntries()>100)
        {
            double mean = hist_time_residual_mpmt[i]->GetMean();
            double sigma = hist_time_residual_mpmt[i]->GetStdDev();
            hist_time_residual_fit->SetBinContent(i,mean);
            hist_time_residual_fit->SetBinError(i,sigma);

            // hist_time_residual_mpmt[i]->Draw("hist");
            // gPad->SaveAs(Form("fig/time_residual_%i.pdf",i));

        }
    }

    fout->cd();

    hist_time_residual_fit->Draw("E");
    gPad->SaveAs(Form("time_residual_fit%s.pdf",suffix));

    hist_end_cap_1->Draw("colz");
    gPad->SaveAs(Form("vertices_1%s.pdf",suffix));
    hist_end_cap_2->Draw("colz");
    gPad->SaveAs(Form("vertices_2%s.pdf",suffix));
    hist_end_cap_3->Draw("colz");
    gPad->SaveAs(Form("vertices_3%s.pdf",suffix));

    hist_bad_hit_ratio->Draw();
    gPad->SaveAs(Form("bad_hit_ratio%s.pdf",suffix));

    gStyle->SetOptStat(0);
    for (int i=0;i<4;i++)
    {
        for (int j=1;j<=hist_QR_Costh[i]->GetNbinsX();j++)
        {
            TH1D* hQR = hist_QR_Costh[i]->ProjectionY("sliceQR", j, j);
            hist_QR_Costh_1D[i]->SetBinContent(j,hQR->GetMean());
            // std::cout<<"slide "<<j<<" nEntries = "<<hQR->GetEntries()<<std::endl;
        }
        hist_QR_Costh_1D[i]->Write();
        hist_QR_Costh_1D[i]->Draw("hist");
        gPad->SaveAs(Form("QR_Costh_1D%s_%i.pdf",suffix,i));
    }
    gStyle->SetOptStat(1);

    for (int i=0;i<106;i++)
    {
        for (int j=1;j<=hist_QR_Costh_mPMT[i]->GetNbinsX();j++)
        {
            TH1D* hY = hist_QR_Costh_mPMT[i]->ProjectionY("sliceY", j, j);
            hist_QR_Costh_mPMT_1D[i]->SetBinContent(j,hY->GetMean());
        }
        hist_QR_Costh_mPMT_1D[i]->Write();
    }

    fout->Close();


    f->Close();
}

void FitAngularResponse()
{
    Read_mPMT_map();

    Fit("/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/mdt_nominal_eff/muon_MC_nominal_selected.root","angular_MC_nominal.root","_MC_nominal");
    Fit("/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/mdt_reduced_eff/muon_MC_selected.root","angular_MC.root","_MC");
    Fit("/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/run1766/muon_1766_selected.root","angular_1766.root","_1766");

    TFile* fNominal = new TFile("angular_MC_nominal.root");
    TFile* fMC = new TFile("angular_MC.root");
    TFile* fData = new TFile("angular_1766.root");

    gStyle->SetOptStat(0);
    std::vector<TH1D*> hNominal(4);
    std::vector<TH1D*> hMC(4);
    std::vector<TH1D*> hData(4);
    for (int i=0;i<4;i++)
    {
        // Plot QR ~ Angular response per mPMT type
        hNominal[i] = (TH1D*)fNominal->Get(Form("hist_QR_Costh_1D_%i",i));
        hMC[i] = (TH1D*)fMC->Get(Form("hist_QR_Costh_1D_%i",i));
        hData[i] = (TH1D*)fData->Get(Form("hist_QR_Costh_1D_%i",i));

        hNominal[i]->SetLineColor(kBlack);
        hMC[i]->SetLineColor(kBlue);
        hData[i]->SetLineColor(kRed);
        hNominal[i]->Draw("hist");
        hMC[i]->Draw("hist same");
        hData[i]->Draw("hist same");
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SaveAs(Form("QR_Costh_All_%i.pdf",i));
    }

    // Plot angular response ratio to nominal MC
    int lineC[4] = {kGreen, kRed, kOrange, kBlue};
    std::vector<std::string> mpmt_type_string{"WUT_exsitu","TRI_exsitu","WUT_insitu","TRI_insitu"};
    // Reduced efficiency MC
    for (int i=0;i<4;i++)
    {
        hMC[i]->Divide(hNominal[i]);
        hMC[i]->SetLineColor(lineC[i]);
        if (i==0) hMC[i]->Draw("hist");
        else hMC[i]->Draw("hist same");
        hMC[i]->GetYaxis()->SetRangeUser(0.4,1);
    }
    gPad->SaveAs("QR_Costh_1D_MC_Ratio.pdf");
    // Data
    for (int i=0;i<4;i++)
    {
        hData[i]->Divide(hNominal[i]);
        hData[i]->SetLineColor(lineC[i]);
        if (i==0) hData[i]->Draw("hist");
        else hData[i]->Draw("hist same");
        hData[i]->GetYaxis()->SetRangeUser(0.2,1.5);
    }
    gPad->SaveAs("QR_Costh_1D_Data_Ratio.pdf");
    // Plot Angular response per mPMT slot ID
    for (int i=0;i<4;i++)
    {
        bool firstData = true;
        double low = 0, high = 1000;
        if (i!=0) high=2000;
        int color_count = 2;
        auto legend = new TLegend(0.1,0.7,0.48,0.9);
        for (int j=0;j<106;j++)
        {
            int mpmt_type_id = mpmt_type[j];
            if (mpmt_type_id==i)
            {
                TH1D* hData = (TH1D*)fData->Get(Form("hist_QR_Costh_1D_Slot_%i",j));
                if (hData->GetMaximum()<10) continue;
                hData->GetYaxis()->SetRangeUser(low,high);
                hData->SetLineColor(kBlack);
                if (hData->GetBinContent(11)<600 && i!=0) 
                {
                    hData->SetLineColor(color_count++);
                    hData->SetLineWidth(3);
                    legend->AddEntry(hData,Form("Slot_%i",j),"l");
                }
                if (firstData) 
                {
                    hData->SetTitle(Form("QR_Costh_1D_Data_%s.pdf",mpmt_type_string[i].c_str()));
                    hData->Draw("hist");
                    firstData = false;
                }
                else hData->Draw("hist same");
            }
        }
        legend->Draw();
        gPad->SaveAs(Form("QR_Costh_1D_Data_%s.pdf",mpmt_type_string[i].c_str()));
    }
}