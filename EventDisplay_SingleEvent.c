R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

void EventDisplay_SingleEvent(const char * fname, int evtID)
{
    gStyle->SetOptStat(0);

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

    std::string prefix = fname;
    if (prefix.find_last_of("/")!=std::string::npos) 
    {   
        prefix = prefix.substr(prefix.find_last_of("/")+1);
    }
    if (prefix.find_last_of("[")!=std::string::npos) 
    {   
        prefix = prefix.substr(0,prefix.find_last_of("["));
    }
    if (prefix.find_last_of("*")!=std::string::npos) 
    {   
        prefix = prefix.substr(0,prefix.find_last_of("*"));
    }
    if (prefix.substr(prefix.length()-5)==".root") 
    {   
        prefix = prefix.substr(0,prefix.length()-5);
    }
    if (prefix.find_last_of("_")!=prefix.length()-1) prefix += ("_");
    std::cout<<"prefix = "<<prefix<<std::endl;

    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    t->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent);

    WCSimRootTrigger* wcsimrootevent;
    // Get vertex and beam direction from the specific event
    t->GetEntry(evtID);
    wcsimrootevent=wcsimrootsuperevent->GetTrigger(0);
    TVector3 vtx(wcsimrootevent->GetVtx(0),wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
    TVector3 BeamDir(((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(0),((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(1),((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(2));
    std::cout<<"BeamDir = "<<BeamDir.x()<<" "<<BeamDir.y()<<" "<<BeamDir.z()<<std::endl;

    // Geometry tree - only need 1 "event"
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
    std::vector<std::vector<double>> pmt_pos(nPMTs_type0);
    std::vector<TVector3> pmt_posT(nPMTs_type0);
    std::vector<std::vector<double>> pmt_dir(nPMTs_type0);
    std::vector<double> pmt_ang(nPMTs_type0);
    std::vector<double> pmt_tof(nPMTs_type0);
    double vg = 2.20027795333758801e8*100/1.e9; // rough speed of light in water in cm/ns
    double max_r = 0, max_z = 0;
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
        pmt_pos[i] = pos;
        pmt_dir[i] = dir;

        TVector3 pmtpos(pos[0],pos[1],pos[2]);
        pmt_posT[i] = pmtpos;
        pmt_ang[i] = BeamDir.Angle(pmtpos-vtx)*180/TMath::Pi();

        pmt_tof[i] = (pmtpos-vtx).Mag()/vg;

        // y-axis is vertical
        if (max_z<fabs(pos[1])) max_z=fabs(pos[1]);
        if (max_r<sqrt(pos[0]*pos[0]+pos[2]*pos[2]))
            if (fabs(pmt.GetOrientation(1))>0.5) max_r = sqrt(pos[0]*pos[0]+pos[2]*pos[2]);
    }

    double barrelCut = max_z-10;
    TH2D* hist_event_display = new TH2D("Charges","Charges",250,-TMath::Pi()*max_r,TMath::Pi()*max_r,250,-max_z-2*max_r,max_z+2*max_r);
    TH2D* hist_event_display_time = new TH2D("TIme","Time",250,-TMath::Pi()*max_r,TMath::Pi()*max_r,250,-max_z-2*max_r,max_z+2*max_r);
    std::vector<std::vector<double>> eventDiplayXY;
    for (int i=0;i<nPMTs_type0;i++)
    {
        // rotation for event display
        double x = -pmt_pos.at(i).at(0);
        double y = pmt_pos.at(i).at(2);
        double z = pmt_pos.at(i).at(1);
        std::vector<double> pmtXY;
        if (fabs(z)<barrelCut) // barrel
        {
            double th = atan2(y,x);
            pmtXY.push_back(-max_r*th);
            pmtXY.push_back(z);
        }
        else if (z>barrelCut) //top
        {
            pmtXY.push_back(-y);
            pmtXY.push_back(max_z+max_r-x);
        }
        else //bot
        {
            pmtXY.push_back(-y);
            pmtXY.push_back(-max_z-max_r+x);
        }
        eventDiplayXY.push_back(pmtXY);
    }

    std::vector<double> pmt_hit(nPMTs_type0,0.);
    std::vector<double> pmt_time(nPMTs_type0,0.);
    TH1D* hist_timetof = new TH1D("DigiTime-TOF","DigiTime-TOF",1000,-20,40);
    TH1D* hist_timetof_true = new TH1D("TrueTime-TOF","TrueTime-TOF",1000,-10,50);
    double time_min = 9999;
    for (long int nev=evtID;nev<evtID+1;nev++)
    {

        std::cout<<"Running "<<nev<<"-th event of total "<<t->GetEntries()<<" events"<<std::endl;

        delete wcsimrootsuperevent;
        wcsimrootsuperevent = 0;  // EXTREMELY IMPORTANT

        t->GetEntry(nev);
        wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

        std::vector<double> triggerInfo = wcsimrootevent->GetTriggerInfo();
        double triggerShift=0, triggerTime=0;
        if(wcsimrootevent->GetTriggerType()!=kTriggerNoTrig && triggerInfo.size()>=3)
        {
            triggerShift = triggerInfo[1];
            triggerTime = triggerInfo[2];
        }

        int nhits = wcsimrootevent->GetNcherenkovdigihits(); 

        // Fill digi hit histogram
        for (int i=0; i< nhits ; i++)
        {
            WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit = (WCSimRootCherenkovDigiHit*) (wcsimrootevent->GetCherenkovDigiHits())->At(i);
            int tubeNumber     = wcsimrootcherenkovdigihit->GetTubeId()-1;
            double peForTube      = wcsimrootcherenkovdigihit->GetQ();
            double time = wcsimrootcherenkovdigihit->GetT()+triggerTime-triggerShift;

            pmt_hit[tubeNumber] += peForTube;
            pmt_time[tubeNumber] = wcsimrootcherenkovdigihit->GetT();
            if (time_min>pmt_time[tubeNumber]) time_min = pmt_time[tubeNumber];

            hist_timetof->Fill(time-pmt_tof[tubeNumber],peForTube);
        }

        // Fill true hit histgram
        int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
        TClonesArray *timeArray = wcsimrootevent->GetCherenkovHitTimes();
        for (int i=0; i< ncherenkovhits ; i++)
        {
            WCSimRootCherenkovHit *wcsimrootcherenkovhit = (WCSimRootCherenkovHit*) (wcsimrootevent->GetCherenkovHits())->At(i);
            int tubeNumber     = wcsimrootcherenkovhit->GetTubeID()-1;
            int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
            int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
            for (int idx = timeArrayIndex; idx<timeArrayIndex+peForTube; idx++)
            {
                WCSimRootCherenkovHitTime * cht = (WCSimRootCherenkovHitTime*) timeArray->At(idx);
                TVector3 endPos(cht->GetPhotonEndPos(0)/10.,cht->GetPhotonEndPos(1)/10.,cht->GetPhotonEndPos(2)/10.); // mm to cm
                TVector3 endDir(cht->GetPhotonEndDir(0),cht->GetPhotonEndDir(1),cht->GetPhotonEndDir(2));
                hist_timetof_true->Fill(cht->GetTruetime()-(endPos-vtx).Mag()/vg);
            }
        }
    }
    
    double time_max = time_min+20;
    for (int i=0;i<nPMTs_type0;i++)
    {
        if (pmt_hit[i]>0)
        {
            hist_event_display->Fill(eventDiplayXY.at(i).at(0),eventDiplayXY.at(i).at(1),pmt_hit[i]);
            if (pmt_time[i]>=time_min && pmt_time[i]<=time_max) hist_event_display_time->Fill(eventDiplayXY.at(i).at(0),eventDiplayXY.at(i).at(1),pmt_time[i]);
            std::cout<<"pmt_time = "<<pmt_time[i]<<std::endl;
        }
    }
    TCanvas* c1 = new TCanvas();

    hist_event_display->Draw("colz");
    double vtx_x = -vtx.x(), vtx_y = vtx.z(), vtx_z = vtx.y();
    // Extrapolate to get muon entrance point
    double entrance_x = vtx_x - BeamDir.x(), entrance_y = vtx_y + BeamDir.z(), entrance_z = vtx_z + BeamDir.y();
    int nSteps_max = 10000;
    int count = 0;
    while (!(sqrt(entrance_x*entrance_x+entrance_y*entrance_y)<max_r && fabs(entrance_z)<max_z)&& count++<nSteps_max)
    {
        entrance_x += -BeamDir.x(); entrance_y += BeamDir.z(); entrance_z += BeamDir.y();
    }
    if (count>=nSteps_max) 
    {
        std::cout<<"Muon not entering detector --> Exit"<<std::endl;
        return;
    }
    // Extrapolate to get muon exit point
    double exit_x = entrance_x - BeamDir.x(), exit_y = entrance_y + BeamDir.z(), exit_z = entrance_z + BeamDir.y();
    while (sqrt(exit_x*exit_x+exit_y*exit_y)<max_r && fabs(exit_z)<max_z)
    {
        exit_x += -BeamDir.x(); exit_y += BeamDir.z(); exit_z += BeamDir.y();
    }
    double evtx, evty;
    if (fabs(entrance_z)<barrelCut)
    {
        double th = atan2(entrance_y,entrance_z);

        evtx = -max_r*th;
        evty = entrance_z;
    }
    else if (entrance_z>barrelCut)
    {
        evtx = -entrance_y;
        evty = max_z+max_r-entrance_x;
    }
    else
    {
        evtx = -entrance_y;
        evty = -max_z-max_r+entrance_x;
    }
    TMarker m1(evtx,evty,29);
    m1.SetMarkerColor(kRed);
    m1.Draw();
    if (fabs(exit_z)<barrelCut)
    {
        double th = atan2(exit_y,exit_x);

        evtx = -max_r*th;
        evty = exit_z;
    }
    else if (exit_z>barrelCut)
    {
        evtx = -exit_y;
        evty = max_z+max_r-exit_x;
    }
    else
    {
        evtx = -exit_y;
        evty = -max_z-max_r+exit_x;
    }
    TMarker m2(evtx,evty,29);
    m2.SetMarkerColor(kBlack);
    m2.Draw();
    c1->SaveAs(Form("%sdisplay_%i.pdf",prefix.c_str(),evtID));

    hist_event_display_time->GetZaxis()->SetRangeUser(time_min,time_max);
    hist_event_display_time->Draw("colz");
    m1.Draw();
    m2.Draw();
    c1->SaveAs(Form("%sdisplay_time_%i.pdf",prefix.c_str(),evtID));

    hist_timetof->GetXaxis()->SetTitle("Digi Time (ns)");
    hist_timetof_true->GetXaxis()->SetTitle("Raw Time (ns)");

    hist_timetof->Draw("hist");
    c1->SaveAs(Form("%stimetof_%i.pdf",prefix.c_str(),evtID));

    hist_timetof_true->Draw("hist");
    //c1->SetLogy();
    c1->SaveAs(Form("%stimetof_true_%i.pdf",prefix.c_str(),evtID));

    f->Close();
    t->Reset();
}