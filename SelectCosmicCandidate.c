R__LOAD_LIBRARY($WCSIM_BUILD_DIR/lib/libWCSimRoot.so)

void SelectCosmicCandidate(int run_num = 1766)
{
    // gStyle->SetOptStat(0);
    
    int n50_thres = 700;

    double max_r = 307.5926/2;
    double max_z = 271.4235/2;

    std::string outfilename = Form("out_%i.root",run_num);
    std::cout<<"outfilename = "<<outfilename<<std::endl;
    std::string offline_files = (Form("/eos/experiment/wcte/data/2025_commissioning/processed_offline_data/production_v0_5/%i/WCTE_offline_R%iS0P*.root",run_num,run_num));
    std::cout<<"offline_files = "<<offline_files<<std::endl;

    std::string geofile = "geofile_NuPRISMBeamTest_16cShort_mPMT.txt";
    ifstream infile(geofile.c_str());
    // skip N header lines
    std::string line;
    for (int i=0; i<5; i++) std::getline(infile, line);

    int cab_id, mpmt_id, pmt_id, cycloc;
    double x, y, z, dx, dy, dz;
    std::vector<std::vector<int>>   mpmt_pmt_cab(106,std::vector<int>(19,-1)); // channel mapping
    std::vector<std::vector<int>>   mpmt_pmt_cyc(106,std::vector<int>(19,-1)); // 0 = top cap, 1 = barrel, 2 = bottom cap
    std::vector<std::vector<TVector3>>   mpmt_pmt_pos(106,std::vector<TVector3>(19));
    std::vector<std::vector<TVector3>>   mpmt_pmt_displayXY(106,std::vector<TVector3>(19)); // PMT position for 2D event display
    while ( infile >> cab_id >> mpmt_id >> pmt_id >> x >> y >> z >> dx >> dy >> dz >> cycloc )
    {
        mpmt_pmt_cab[mpmt_id][pmt_id-1] = cab_id;
        mpmt_pmt_cyc[mpmt_id][pmt_id-1] = cycloc;
        mpmt_pmt_pos[mpmt_id][pmt_id-1] = TVector3(x,y,z);
        // std::cout<<x<<" "<<y<<" "<<z<<"\n";
        y = -mpmt_pmt_pos[mpmt_id][pmt_id-1].x();
        x =  mpmt_pmt_pos[mpmt_id][pmt_id-1].z();
        z =  mpmt_pmt_pos[mpmt_id][pmt_id-1].y();
        if (cycloc==1)
        {
            double th = atan2(y,x);
            mpmt_pmt_displayXY[mpmt_id][pmt_id-1] = TVector3(-max_r*th,z,0);
        }
        else if (cycloc==0)
        {
            mpmt_pmt_displayXY[mpmt_id][pmt_id-1] = TVector3(-y,max_z+max_r-x,0);
        }
        else
        {
            mpmt_pmt_displayXY[mpmt_id][pmt_id-1] = TVector3(-y,-max_z-max_r+x,0);
        }
    }

    TFile * outfile = new TFile(outfilename.c_str(),"RECREATE");
    TTree* fWCSimT = new TTree("wcsimT", "Mockup WCSim Tree for fiTQun");
    WCSimRootEvent* fSpEvt = new WCSimRootEvent();
    fSpEvt->Initialize();
    fWCSimT->Branch("wcsimrootevent", "WCSimRootEvent", &fSpEvt);
    // TH2D* hist_event_display = new TH2D("Charges","Charges",250,-TMath::Pi()*max_r,TMath::Pi()*max_r,250,-max_z-2*max_r,max_z+2*max_r);

    std::cout<<"Load WCTE data: "<<offline_files<<std::endl;
    TChain *chain = new TChain("WCTEReadoutWindows");
    chain->Add(offline_files.c_str());
    long int nevent = ((int)chain->GetEntries());
    printf("nevent %ld\n",nevent);
    std::vector<float> * hit_pmt_charges = 0;
    std::vector<double> * hit_pmt_times = 0;
    std::vector<int> * hit_mpmt_slot_ids = 0;
    std::vector<int> * hit_pmt_position_ids = 0;
    std::vector<int> * hit_mpmt_card_ids = 0;
    // std::vector<bool> * hit_pmt_has_time_constant = 0;
    int window_data_quality;
    std::vector<int> * hit_pmt_readout_mask = 0;
    chain->SetBranchAddress("hit_pmt_charges",&hit_pmt_charges);
    chain->SetBranchAddress("hit_pmt_calibrated_times",&hit_pmt_times);
    chain->SetBranchAddress("hit_mpmt_slot_ids",&hit_mpmt_slot_ids);
    chain->SetBranchAddress("hit_pmt_position_ids",&hit_pmt_position_ids);
    chain->SetBranchAddress("hit_mpmt_card_ids",&hit_mpmt_card_ids);
    // chain->SetBranchAddress("hit_pmt_has_time_constant",&hit_pmt_has_time_constant);
    chain->SetBranchAddress("window_data_quality",&window_data_quality);
    chain->SetBranchAddress("hit_pmt_readout_mask",&hit_pmt_readout_mask);

    int cosmics_count = 0;
    int count1pc = nevent/100;
    if (count1pc==0) count1pc=1;
    for (long int i=0;i<nevent;i++)
    {
        if (i%count1pc==0) std::cout<<"Running "<<i<<"-th event of total "<<nevent<<" events"<<std::endl;
        chain->GetEntry(i);
        // printf("first hit %f %f\n",(*hit_pmt_charges)[0],(*hit_pmt_times)[0]);

        // 0= no potential problems found, 1= trigger effected by 67ms issue, 2= slow control data quality issues identified (e.g. dropped packets)
        if (window_data_quality!=0) continue;

        int len = (*hit_mpmt_card_ids).size();

        // Select valid hits and create an index array
        std::vector<double> times;
        std::vector<size_t> indices;
        for (int j=0;j<len;j++)
        {
            indices.push_back(j);
            times.push_back(-9999);

            if ((*hit_mpmt_card_ids)[j]>=120) continue;
            if ((*hit_pmt_charges)[j]>10000) continue;
            // if ((*hit_pmt_has_time_constant)[j]==0) continue;
            if ((*hit_pmt_readout_mask)[j]!=0) continue;

            int mpmtid = (*hit_mpmt_slot_ids)[j];
            int  pmtid = (*hit_pmt_position_ids)[j];
            times[j] = (*hit_pmt_times)[j];
        }

        // Sort the index array using times as a reference
        std::sort(indices.begin(), indices.end(),
                [&times](size_t i1, size_t i2) { return times[i1] < times[i2]; });
        // Sort times itself
        std::sort(times.begin(),times.end());

        std::vector<std::vector<double>> trigger_times;
        for (int j=0;j<len;j++)
        {
            double t = times[j];
            if (t<=-9999) continue;
            // lower_bound returns an iterator to the first element not less than value
            auto it = std::lower_bound(times.begin(), times.end(), t+50); // n50 search
            int count = it - times.begin() -j;
            if (count > n50_thres) // select nhit cluster 
            {
                // further include hits that are less than 50 ns apart
                int last_hit_pos = it - times.begin() - 1;
                for (int k=last_hit_pos+1;k<len;k++)
                {
                    if (times[k]-times[last_hit_pos]<50) last_hit_pos = k;
                    else break;
                }

                std::vector<double> pmt_hit(3,0.); // hit count in three (top,barrel,bottom) sectors
                double totQ = 0;
                // hist_event_display->Reset();
                // Prepare WCSim event structure
                fSpEvt->ReInitialize();
                // one trigger per cosmic candidate
                WCSimRootTrigger* anEvent = fSpEvt->GetTrigger(0);
                TriggerType_t trigType = kTriggerNoTrig;
                float hitTimeOffset = 0;
                vector<Double_t> info(1, last_hit_pos-j+1); // nhit
                info.push_back(hitTimeOffset);
                info.push_back(times[j]); // trigger time
                anEvent->SetTriggerInfo(trigType, info);
                for (int k=j;k<=last_hit_pos;k++)
                {
                    int idx = indices[k];
                    mpmt_id = (*hit_mpmt_slot_ids)[idx];
                    pmt_id = (*hit_pmt_position_ids)[idx];
                    cycloc = mpmt_pmt_cyc[mpmt_id][pmt_id];

                    pmt_hit[cycloc] += (*hit_pmt_charges)[idx];
                    totQ += (*hit_pmt_charges)[idx];

                    double hitq = (*hit_pmt_charges)[idx]/135.;
                    double hitt = times[k]-times[j];
                    int tubeID = mpmt_pmt_cab[mpmt_id][pmt_id];
                    std::vector<int> true_pe_comp;
                    // Fill hit in trigger
                    anEvent->AddCherenkovDigiHit(hitq,
                                                hitt,
                                                tubeID,
                                                mpmt_id,
                                                pmt_id+1,
                                                true_pe_comp);

                    // hist_event_display->Fill(mpmt_pmt_displayXY[mpmt_id][pmt_id].x(),mpmt_pmt_displayXY[mpmt_id][pmt_id].y(),(*hit_pmt_charges)[idx]);
                }
                anEvent->SetNumDigitizedTubes(last_hit_pos-j+1);
                anEvent->SetSumQ(totQ);
                anEvent->GetHeader()->SetDate(int(times[j]));
                j = last_hit_pos;
                
                // hit ratio cuts, based on MC
                if (pmt_hit[0]/totQ>0.07) continue;
                if (pmt_hit[1]/totQ<0.45 || pmt_hit[1]/totQ>0.7) continue;
                if (pmt_hit[2]/totQ<0.25 || pmt_hit[2]/totQ>0.5) continue;
                fWCSimT->Fill();


                // hist_event_display->Draw("colz");
                // gPad->SetLogz();
                // gPad->SaveAs(Form("fig/cosmics_%i.pdf",cosmics_count++));
                // std::cout<<"cosmics candidate!!!\n";
            }
        }
    }

    // Copy metadata
    const vector<string> listWCRootCopyTree{"wcsimGeoT","Settings","wcsimRootOptionsT"};
    TFile* wcsimFile = new TFile("wcsim_dummy.root");
    for (auto s : listWCRootCopyTree)
    {
        TTree* tin = (TTree*)wcsimFile->Get(s.c_str());
        outfile->cd();
        TTree *tout = tin->CloneTree(-1, "fast");
        tout->Write();
    }
    wcsimFile->Close();

    outfile->cd();
    fWCSimT->Write();
    outfile->Close();
}