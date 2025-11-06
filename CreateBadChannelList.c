#include "json.hpp"
using json = nlohmann::json;

void CreateBadChannelList(int run_num=1766)
{
    // first default to mask all PMTs
    std::vector<bool> pmt_mask_list(1843,true);

    // Read channel map between WCSim ID and (mPMT_slot_id,pmt_position_id)
    std::string geofile = "geofile_NuPRISMBeamTest_16cShort_mPMT.txt";
    ifstream infile(geofile.c_str());
    // skip N header lines
    std::string line;
    for (int i=0; i<5; i++) std::getline(infile, line);

    int cab_id, mpmt_id, pmt_id, cycloc;
    double x, y, z, dx, dy, dz;
    std::vector<std::vector<int>>   mpmt_pmt_cab(106,std::vector<int>(19,-1));
    while ( infile >> cab_id >> mpmt_id >> pmt_id >> x >> y >> z >> dx >> dy >> dz >> cycloc )
    {
        // mpmt_id is from 0 - 105, while pmt_id is from 1-19
        // cab_id is starting from 1
        mpmt_pmt_cab[mpmt_id][pmt_id-1] = cab_id-1;
    }

    // Un-mask good channel list from metadata
    std::string meta_data_file = (Form("/eos/experiment/wcte/data/2025_commissioning/processed_offline_data/production_v0_5/%i/run_%i_meta_data_json.json",run_num,run_num));
    std::cout<<"meta_data_file = "<<meta_data_file<<std::endl;

    std::ifstream ifs(meta_data_file.c_str());
    json data = json::parse(ifs);

    for (int i=0;i<data.at("good_wcte_pmts").size();i++)
    {
        // channel_id = 100*mPMT_slot_id + pmt_position_id
        int channel_id = data.at("good_wcte_pmts").at(i);
        mpmt_id = channel_id/100;
        pmt_id = channel_id%100;
        // std::cout<<"Un-mask "<<mpmt_id<<" "<<pmt_id<<" "<<mpmt_pmt_cab[mpmt_id][pmt_id]<<std::endl;
        pmt_mask_list[mpmt_pmt_cab[mpmt_id][pmt_id]]=false;
    }

    // Further, Check run data and flag channels with no hits
    std::string offline_files = (Form("/eos/experiment/wcte/data/2025_commissioning/processed_offline_data/production_v0_5/%i/WCTE_offline_R%iS0P*.root",run_num,run_num));
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
    std::vector<bool> active_pmts(1843,false);
    // nevent = 10000;
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

        // Select valid hits 
        for (int j=0;j<len;j++)
        {
            if ((*hit_mpmt_card_ids)[j]>=120) continue;
            if ((*hit_pmt_charges)[j]>10000) continue;
            // if ((*hit_pmt_has_time_constant)[j]==0) continue;
            if ((*hit_pmt_readout_mask)[j]!=0) continue;

            int mpmtid = (*hit_mpmt_slot_ids)[j];
            int  pmtid = (*hit_pmt_position_ids)[j];
            active_pmts[mpmt_pmt_cab[mpmtid][pmtid]]=true;
        }
    }

    std::cout<<"Create mask file: "<<Form("maskFile_%i.txt",run_num)<<std::endl;
    std::ofstream maksFile;   // File for text output
    maksFile.open(Form("maskFile_%i.txt",run_num), std::ios::out);
    for (int i=0;i<1843;i++)
    {
        // mask PMT with no hits
        if (!pmt_mask_list[i]&&!active_pmts[i]) pmt_mask_list[i]=true;

        if (pmt_mask_list[i]) maksFile << i << std::endl;
    }
    maksFile.close();
}