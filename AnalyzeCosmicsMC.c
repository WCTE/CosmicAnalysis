R__LOAD_LIBRARY($WCSIM_BUILD_DIR/lib/libWCSimRoot.so)

void AnalyzeCosmicsMC(
    const char * fname= "/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/mdt/mdt_cosmics_00*.root", 
    const char * fqname="/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/mdt/fiTQun/fq_00*.root", 
    bool topdownonly = false, bool applyselection=true, const char * maks_file_name="maskFile_1766.txt"
)
{
    // gStyle->SetOptStat(0);

    TChain *t = new TChain("wcsimT");
    t->Add(fname);
    std::string single_file_name = t->GetFile()->GetName();
    // Get the first file for geometry
    TFile *f = TFile::Open(single_file_name.c_str());
    if (!f->IsOpen()){
        std::cout << "Error, could not open input file: " << single_file_name << std::endl;
        return -1;
    }
    if (!f->IsOpen()) return;

    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    t->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent);

    TChain *tFQ = new TChain("fiTQun");
    tFQ->Add(fqname);
    t->AddFriend(tFQ);
    int fqnse;
    float fq1rt0[100][7];
    float fq1rpos[100][7][3];
    float fq1rdir[100][7][3];
    float fq1rnll[100][7];
    float fq1rmom[100][7];
    t->SetBranchAddress("fqnse",&fqnse);
    t->SetBranchAddress("fq1rt0",&fq1rt0);
    t->SetBranchAddress("fq1rpos",&fq1rpos);
    t->SetBranchAddress("fq1rdir",&fq1rdir);
    t->SetBranchAddress("fq1rnll",&fq1rnll);
    t->SetBranchAddress("fq1rmom",&fq1rmom);

    WCSimRootTrigger* wcsimrootevent;

    // Geometry tree - only need 1 "event"
    double max_r = 307.5926/2;
    double max_z = 271.4235/2;
    WCSimRootGeom *geo = 0;
    TTree *geotree = (TTree*)f->Get("wcsimGeoT");
    geotree->SetBranchAddress("wcsimrootgeom", &geo);
    std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
    if (geotree->GetEntries() == 0) {
        exit(9);
    }
    geotree->GetEntry(0);
    int nPMTs_type0=geo->GetWCNumPMT();
    std::cout << "geo has " << nPMTs_type0 << " PMTs" << std::endl;
    std::vector<TVector3> pmt_posT(nPMTs_type0);
    std::vector<int> pmt_CylLoc(nPMTs_type0);
    std::vector<int> mPMTNo(nPMTs_type0);
    std::vector<int> mPMT_PMTNo(nPMTs_type0);
    for (int i=0;i<nPMTs_type0;i++) 
    {
        WCSimRootPMT pmt;
        pmt = geo->GetPMT(i);
        std::vector<double> pos(3);
        std::vector<double> dir(3);
        for(int j=0;j<3;j++){
            pos[j] = pmt.GetPosition(j);
            dir[j] = pmt.GetOrientation(j);
        }
        mPMTNo[i] = pmt.GetmPMTNo();
        mPMT_PMTNo[i] = pmt.GetmPMT_PMTNo();
        pmt_CylLoc[i] = pmt.GetCylLoc();

        // std::cout<<i<<" "<<mPMTNo[i]<<" "<<mPMT_PMTNo[i]<<" "<<pmt_CylLoc[i]<<"\n";

        TVector3 pmtpos(pos[0],pos[1],pos[2]);
        pmt_posT[i] = pmtpos;
    }

    ifstream maskFile(maks_file_name);
    int cab_id;
    std::vector<bool> mask_list(nPMTs_type0,false);
    while (maskFile >> cab_id)
    {
        mask_list[cab_id] = true;
        // std::cout<<"cab_id "<<cab_id<<" is masked\n";
    }


    TH1D* hist_top = new TH1D("hist_top","hist_top",100,0,1);
    TH1D* hist_barrel = new TH1D("hist_barrel","hist_barrel",100,0,1);
    TH1D* hist_bottom = new TH1D("hist_bottom","hist_bottom",100,0,1);
    TH1D* hist_totQ = new TH1D("hist_totQ","hist_totQ",100,0,10000);
    TH1D* hist_nhit = new TH1D("hist_nhit","hist_nhit",100,0,1843);
    TH1D* hist_pe = new TH1D("hist_pe","hist_pe",1000,0,30);
    TH1D* hist_dir = new TH1D("hist_dir","hist_dir",100,-1,1);

    TH2D* hist_top_vs_dir = new TH2D("hist_top_vs_dir","hist_top_vs_dir",100,0,1,100,0,1);
    TH2D* hist_barrel_vs_dir = new TH2D("hist_barrel_vs_dir","hist_barrel_vs_dir",100,0,1,100,0,1);
    TH2D* hist_bottom_vs_dir = new TH2D("hist_bottom_vs_dir","hist_bottom_vs_dir",100,0,1,100,0,1);
    TH2D* hist_totQ_vs_dir = new TH2D("hist_totQ_vs_dir","hist_totQ_vs_dir",100,0,1,100,0,20000);
    TH2D* hist_nhit_vs_dir = new TH2D("hist_nhit_vs_dir","hist_nhit_vs_dir",100,0,1,100,0,1843);

    TH2D* hist_entrance = new TH2D("Entrance_Pos","Entrance_Pos",100,-TMath::Pi()*max_r,TMath::Pi()*max_r,100,-max_z-2*max_r,max_z+2*max_r);
    TH2D* hist_exit = new TH2D("Exit_Pos","Exit_Pos",100,-TMath::Pi()*max_r,TMath::Pi()*max_r,100,-max_z-2*max_r,max_z+2*max_r);
    TH1D* hist_dist = new TH1D("Dist","Dist",100,0,500);

    TH2D* hist_fqentrance = new TH2D("FQ_Entrance_Pos","FQ_Entrance_Pos",100,-TMath::Pi()*max_r,TMath::Pi()*max_r,100,-max_z-2*max_r,max_z+2*max_r);
    TH2D* hist_fqexit = new TH2D("FQ_Exit_Pos","FQ_Exit_Pos",100,-TMath::Pi()*max_r,TMath::Pi()*max_r,100,-max_z-2*max_r,max_z+2*max_r);
    TH2D* hist_fqdirxy = new TH2D("FQ_Direction_Y_vs_X","FQ_Direction_Y_vs_X",100,-1,1,100,-1,1);
    TH1D* hist_fqdist = new TH1D("FQ_Dist","FQ_Dist",100,0,500);
    TH1D* hist_fqdirz = new TH1D("FQ_Dirz","FQ_Dirz",100,0.5,1);
    TH1D* hist_fqposz = new TH1D("pos_z","pos_z",100,-200,200);
    TH1D* hist_fqmom = new TH1D("FQ_Mom","FQ_Mom",100,0,1000);
    TH2D* hist_fqmom_vs_Q = new TH2D("FQ_Mom_vs_Q","FQ_Mom_vs_Q",100,0,10000,100,400,1000);
    TH1D* hist_fq_charge_over_dist = new TH1D("FQ_Charge_Over_Dist","FQ_Charge_Over_Dist",100,0,50);
    TH1D* hist_fqnll = new TH1D("FQ_mu_e_likelihood","FQ_mu_e_likelihood",100,-1000,1000);
    TH1D* hist_fq_hit_time = new TH1D("FQ_Hit_Time","FQ_Hit_Time",100,-10,40);

    TH1D* hist_fqentrance_diff = new TH1D("FQ_Entrance_Pos_Diff","FQ_Entrance_Pos_Diff",100,0,100);
    TH1D* hist_fqexit_diff = new TH1D("FQ_Exit_Pos_Diff","FQ_Exit_Pos_Diff",100,0,100);
    TH1D* hist_fqdir_diff = new TH1D("FQ_Dir_Diff","FQ_Dir_Diff",100,-1,1);
    TH1D* hist_fqdist_diff = new TH1D("FQ_Dist_Diff","FQ_Dist_Diff",100,-100,100);

    int count1pc = t->GetEntries()/100;
    if (count1pc==0) count1pc=1;
    int count_selected = 0;
    int count_topdown = 0;
    for (long int nev=0;nev<t->GetEntries();nev++)
    {
        if (nev%(count1pc)==0) std::cout<<"Running "<<nev<<"-th event of total "<<t->GetEntries()<<" events"<<std::endl;

        delete wcsimrootsuperevent;
        wcsimrootsuperevent = 0;  // EXTREMELY IMPORTANT

        t->GetEntry(nev);
        wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

        bool topdown = false;
        double dirz = -999;
        double displayTopX,displayTopY,displayBotX,displayBotY;
        TVector3 entrance_truth, exit_truth, dir_truth;
        for (int i=0;i<wcsimrootevent->GetNtrack();i++)
        {
            WCSimRootTrack* trk = (WCSimRootTrack*)wcsimrootevent->GetTracks()->At(i);
            if (trk->GetId()==1) // select primary track to get vertex and direction
            {
                dirz = -trk->GetDir(1);

                std::vector<std::vector<float>> bp = trk->GetBoundaryPoints();
                std::vector<int> bt = trk->GetBoundaryTypes();
                int last_idx = bp.size()-1;
                if (last_idx<0) break; // muon does not enter the ID
                bool topcap = false;
                TVector3 topPt(bp[0][2]/10.,-bp[0][0]/10.,bp[0][1]/10.);
                displayTopX = -topPt.y();
                displayTopY = max_z+max_r-topPt.x();
                if (bt[0]==1 && topPt.z()>max_z) topcap=true; // enter through top cap blacksheet
                else if (bt[0]==2 && topPt.z()>max_z-10) topcap=true; // enter through top cap mPMT
                else
                {
                    double th = atan2(topPt.y(),topPt.x());
                    displayTopX = -max_r*th;
                    displayTopY = topPt.z();
                }
                entrance_truth = topPt;
                TVector3 botPt(bp[last_idx][2]/10.,-bp[last_idx][0]/10.,bp[last_idx][1]/10.);
                displayBotX = -botPt.y();
                displayBotY = -max_z-max_r+botPt.x();
                bool botcap = false;
                if (bt[last_idx]==1 && botPt.z()<-max_z) botcap=true; // exit through bottom cap blacksheet
                else if (bt[last_idx]==2 && botPt.z()<-max_z+10) botcap=true; // exit through bottom cap mPMT
                else
                {
                    double th = atan2(botPt.y(),botPt.x());
                    displayBotX = -max_r*th;
                    displayBotY = botPt.z();
                }
                exit_truth = botPt;

                dir_truth = (exit_truth-entrance_truth).Unit();

                if (topcap && botcap) topdown=true;

                // if (topdown) std::cout<<bp[0][0]<<" "<<bp[0][1]<<" "<<bp[0][2]<<" "<<bp[last_idx][0]<<" "<<bp[last_idx][1]<<" "<<bp[last_idx][2]<<"\n";
                break;
            }
        }

        // skip bad events
        if ((entrance_truth-exit_truth).Mag()<50)  continue;

        // select true top-down events if necessary
        if (topdownonly && !topdown) continue;

        std::vector<double> triggerInfo = wcsimrootevent->GetTriggerInfo();
        double triggerShift=0, triggerTime=0;
        if(wcsimrootevent->GetTriggerType()!=kTriggerNoTrig && triggerInfo.size()>=3)
        {
            triggerShift = triggerInfo[1];
            triggerTime = triggerInfo[2];
        }

        TVector3 fqVtx(fq1rpos[0][2][2],-fq1rpos[0][2][0],fq1rpos[0][2][1]);
        TVector3 fqDir(fq1rdir[0][2][2],-fq1rdir[0][2][0],fq1rdir[0][2][1]);
        bool fqtopdown = true;
        double topcapdist = (max_z-fqVtx.z())/fqDir.z();
        TVector3 topPt = fqVtx+topcapdist*fqDir;
        double displayTopXFQ, displayTopYFQ;
        if (topPt.Perp()>max_r)
        {
            while (topPt.Perp()>max_r && fabs(topPt.z())<=max_z) topPt = topPt+fqDir;
            double th = atan2(topPt.y(),topPt.x());
            displayTopXFQ = -max_r*th;
            displayTopYFQ = topPt.z();
            fqtopdown = false;
        }
        else
        {
            displayTopXFQ = -topPt.y();
            displayTopYFQ = max_z+max_r-topPt.x();
        }
        double botcapdist = (-max_z-fqVtx.z())/fqDir.z();
        TVector3 botPt = fqVtx+botcapdist*fqDir;
        double displayBotXFQ, displayBotYFQ;
        if (botPt.Perp()>max_r)
        {
            while (botPt.Perp()>max_r && fabs(botPt.z())<=max_z) botPt = botPt-fqDir;
            double th = atan2(botPt.y(),botPt.x());
            displayBotXFQ = -max_r*th;
            displayBotYFQ = botPt.z();
            fqtopdown = false;
        }
        else
        {
            displayBotXFQ = -botPt.y();
            displayBotYFQ = -max_z-max_r+botPt.x();
        }

        int nhits = 0;
        double totQ = 0;
        std::vector<double> pmt_hit(3,0.);
        std::vector<double> tof_corrected_time;
        std::vector<double> pmt_charge;
        for (int i=0; i< wcsimrootevent->GetNcherenkovdigihits() ; i++)
        {
            WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit = (WCSimRootCherenkovDigiHit*) (wcsimrootevent->GetCherenkovDigiHits())->At(i);
            int tubeNumber     = wcsimrootcherenkovdigihit->GetTubeId()-1;
            double peForTube      = wcsimrootcherenkovdigihit->GetQ();
            double time = wcsimrootcherenkovdigihit->GetT()+triggerTime-triggerShift; // conversion from local digi time to global time

            if (mask_list[tubeNumber]) continue;
            if (time>50) continue;
            //if (peForTube>30) continue;

            hist_pe->Fill(peForTube);

            int cycloc = pmt_CylLoc[tubeNumber];
            pmt_hit[cycloc] += peForTube;
            nhits++;
            totQ += peForTube;

            double vg = 2.20027795333758801e8*100/1.e9; // rough speed of light in water in cm/ns
            double vg_mu = 30; // rough speed of light in water in cm/ns
            // Residual time
            time -= fq1rt0[0][2] + (max_z-fqVtx.z())/fqDir.z()/vg_mu // fitted entrance time
                    + (pmt_posT[tubeNumber]-TVector3(fq1rpos[0][2][0],fq1rpos[0][2][1],fq1rpos[0][2][2])).Mag()/vg; // photon time of flight from entrance point
            
            tof_corrected_time.push_back(time);
            pmt_charge.push_back(peForTube);
        }

        // selection based on observables
        if (applyselection)
        {
            if (pmt_hit[0]/totQ>0.07) continue;
            if (pmt_hit[1]/totQ<0.45 || pmt_hit[1]/totQ>0.7) continue;
            if (pmt_hit[2]/totQ<0.25 || pmt_hit[2]/totQ>0.5) continue;
            //if (totQ<6000 || totQ>10000) continue;
            if (nhits<700) continue;
            if (fqDir.z()>0) continue;
            if (!fqtopdown) continue;
        }
        hist_nhit->Fill(nhits); hist_nhit_vs_dir->Fill(dirz,nhits);
        hist_top->Fill(pmt_hit[0]/totQ); hist_top_vs_dir->Fill(dirz,pmt_hit[0]/totQ);
        hist_barrel->Fill(pmt_hit[1]/totQ);  hist_barrel_vs_dir->Fill(dirz,pmt_hit[1]/totQ);
        hist_bottom->Fill(pmt_hit[2]/totQ);  hist_bottom_vs_dir->Fill(dirz,pmt_hit[2]/totQ);
        hist_totQ->Fill(totQ); hist_totQ_vs_dir->Fill(dirz,totQ);
        hist_dir->Fill(dirz);
        hist_fqposz->Fill(fqVtx.z());
        hist_entrance->Fill(displayTopX,displayTopY);
        hist_exit->Fill(displayBotX,displayBotY);
        hist_dist->Fill((entrance_truth-exit_truth).Mag()); 
        hist_fqmom->Fill(fq1rmom[0][2]);
        hist_fqmom_vs_Q->Fill(totQ,fq1rmom[0][2]);
        hist_fqnll->Fill(fq1rnll[0][1]-fq1rnll[0][2]);

        hist_fqdirxy->Fill(fqDir.x(),fqDir.y());
        hist_fqdirz->Fill(-fqDir.z());
        hist_fqentrance->Fill(displayTopXFQ,displayTopYFQ);
        hist_fqexit->Fill(displayBotXFQ,displayBotYFQ);
        hist_fqdist->Fill((topPt-botPt).Mag());

        hist_fqentrance_diff->Fill((topPt-entrance_truth).Mag());
        hist_fqexit_diff->Fill((botPt-exit_truth).Mag());
        hist_fqdir_diff->Fill(fqDir.Dot(dir_truth));
        hist_fqdist_diff->Fill((topPt-botPt).Mag()-(entrance_truth-exit_truth).Mag());
        hist_fq_charge_over_dist->Fill(totQ/(topPt-botPt).Mag());
        for (int j=0;j<tof_corrected_time.size();j++) hist_fq_hit_time->Fill(tof_corrected_time[j],pmt_charge[j]);

        count_selected++;
        if (topdown) count_topdown++;
    }

    std::cout<<"Number of selected muons = "<<count_selected<<std::endl;
    std::cout<<"Number of true top-down muons = "<<count_topdown<<std::endl;

    std::string suffix_string = "";
    if (topdownonly) suffix_string += "_topdown";
    if (applyselection) suffix_string += "_selected";
    char* suffix = (char*)suffix_string.c_str();
    if (gSystem->AccessPathName("fig")) gSystem->mkdir("fig");

    hist_top->Draw();
    gPad->SaveAs(Form("fig/pmt_top%s.pdf",suffix));
    hist_barrel->Draw();
    gPad->SaveAs(Form("fig/pmt_barrel%s.pdf",suffix));
    hist_bottom->Draw();
    gPad->SaveAs(Form("fig/pmt_bottom%s.pdf",suffix));
    hist_totQ->Draw();
    gPad->SaveAs(Form("fig/pmt_totQ%s.pdf",suffix));
    hist_nhit->Draw();
    gPad->SaveAs(Form("fig/pmt_nhit%s.pdf",suffix));
    hist_pe->Draw();
    gPad->SaveAs(Form("fig/pmt_pe%s.pdf",suffix));
    // cumulative histogram
    TH1D *hqcum = (TH1D*)hist_pe->Clone("hist_pe_cum");
    hqcum->Reset(); // keep binning but empty content
    hqcum->SetTitle("PE Cumulative distribution");

    double totalq = hist_pe->Integral(0, hist_pe->GetNbinsX()+1); // include under/overflow
    double sumq = 0.0;

    for (int i = 1; i <= hist_pe->GetNbinsX(); ++i) {
        sumq += hist_pe->GetBinContent(i);
        hqcum->SetBinContent(i, sumq / totalq); // normalized CDF
    }
    hqcum->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SaveAs(Form("fig/pmt_pe_cum%s.pdf",suffix));

    hist_top_vs_dir->Draw("colz");
    gPad->SaveAs(Form("fig/pmt_top_vs_dir%s.pdf",suffix));
    hist_barrel_vs_dir->Draw("colz");
    gPad->SaveAs(Form("fig/pmt_barrel_vs_dir%s.pdf",suffix));
    hist_bottom_vs_dir->Draw("colz");
    gPad->SaveAs(Form("fig/pmt_bottom_vs_dir%s.pdf",suffix));
    hist_totQ_vs_dir->Draw("colz");
    gPad->SaveAs(Form("fig/pmt_totQ_vs_dir%s.pdf",suffix));
    hist_nhit_vs_dir->Draw("colz");
    gPad->SaveAs(Form("fig/pmt_nhit_vs_dir%s.pdf",suffix));

    hist_dir->Draw();
    gPad->SaveAs(Form("fig/muon_dir%s.pdf",suffix));
    hist_entrance->Draw("colz");
    gPad->SaveAs(Form("fig/muon_entrance%s.pdf",suffix));
    hist_exit->Draw("colz");
    gPad->SaveAs(Form("fig/muon_exit%s.pdf",suffix));
    hist_dist->Draw();
    gPad->SaveAs(Form("fig/muon_dist%s.pdf",suffix));

    // cumulative histogram
    TH1D *hcum = (TH1D*)hist_dir->Clone("hist_dir_cum");
    hcum->Reset(); // keep binning but empty content
    hcum->SetTitle("Dirz Cumulative distribution");

    double total = hist_dir->Integral(0, hist_dir->GetNbinsX()+1); // include under/overflow
    double sum = 0.0;

    for (int i = 1; i <= hist_dir->GetNbinsX(); ++i) {
        sum += hist_dir->GetBinContent(i);
        hcum->SetBinContent(i, sum / total); // normalized CDF
    }
    hcum->Draw();
    gPad->SaveAs(Form("fig/muon_dir_cum%s.pdf",suffix));


    hist_fqentrance->Draw("colz");
    gPad->SaveAs(Form("fig/fq_entrance_pos%s.pdf",suffix));
    hist_fqexit->Draw("colz");
    gPad->SaveAs(Form("fig/fq_exit_pos%s.pdf",suffix));
    hist_fqdirxy->Draw("colz");
    gPad->SaveAs(Form("fig/fq_dirxy%s.pdf",suffix));
    hist_fqdist->Draw("");
    gPad->SaveAs(Form("fig/fq_dist%s.pdf",suffix));
    hist_fqmom->Draw("");
    gPad->SaveAs(Form("fig/fq_mom%s.pdf",suffix));

    hist_fqentrance_diff->Draw();
    gPad->SaveAs(Form("fig/fq_entrance_diff%s.pdf",suffix));
    hist_fqexit_diff->Draw();
    gPad->SaveAs(Form("fig/fq_exit_diff%s.pdf",suffix));
    hist_fqdir_diff->Draw();
    gPad->SaveAs(Form("fig/fq_dir_diff%s.pdf",suffix));
    hist_fqdist_diff->Draw();
    gPad->SaveAs(Form("fig/fq_dist_diff%s.pdf",suffix));
    hist_fqmom_vs_Q->Draw("colz");
    gPad->SaveAs(Form("fig/fq_mom_vs_Q%s.pdf",suffix));
    hist_fq_charge_over_dist->Draw();
    gPad->SaveAs(Form("fig/fq_charge_over_dist%s.pdf",suffix));
    hist_fqdirz->Draw();
    gPad->SaveAs(Form("fig/fq_dirz%s.pdf",suffix));
    hist_fqposz->Draw();
    gPad->SaveAs(Form("fig/fq_posz%s.pdf",suffix));
    hist_fqnll->Draw();
    gPad->SaveAs(Form("fig/fq_nll%s.pdf",suffix));
    hist_fq_hit_time->Draw("hist");
    gPad->SaveAs(Form("fig/fq_tof%s.pdf",suffix));
}