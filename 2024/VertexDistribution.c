R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

void VertexDistribution(const char * fname)
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

        // y-axis is vertical
        if (max_z<fabs(pos[1])) max_z=fabs(pos[1]);
        if (max_r<sqrt(pos[0]*pos[0]+pos[2]*pos[2]))
            if (fabs(pmt.GetOrientation(1))>0.5) max_r = sqrt(pos[0]*pos[0]+pos[2]*pos[2]);
    }

    int nBins = 150;
    double axisScale = 1.1;
    int nSteps_max = 10000;

    TH3D* hist_vertices = new TH3D("Vertices","Vertices",nBins,-max_r*axisScale,max_r*axisScale,nBins,-max_r*axisScale,max_r*axisScale,nBins,-max_z*axisScale,max_z*axisScale);
    TH2D* hist_vertices_xy = new TH2D("VerticesXY","VerticesXY",nBins,-max_r*axisScale,max_r*axisScale,nBins,-max_r*axisScale,max_r*axisScale);
    TH2D* hist_vertices_yz = new TH2D("VerticesYZ","VerticesYZ",nBins,-max_r*axisScale,max_r*axisScale,nBins,-max_z*axisScale,max_z*axisScale);
    TH2D* hist_vertices_zx = new TH2D("VerticesZX","VerticesZX",nBins,-max_z*axisScale,max_z*axisScale,nBins,-max_r*axisScale,max_r*axisScale);
    for (long int nev=0;nev<t->GetEntries();nev++)
    {
        if (nev%(t->GetEntries()/100)==0) std::cout<<"Running "<<nev<<"-th event of total "<<t->GetEntries()<<" events"<<std::endl;

        delete wcsimrootsuperevent;
        wcsimrootsuperevent = 0;  // EXTREMELY IMPORTANT

        t->GetEntry(nev);
        wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
        TVector3 BeamDir(((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(0),((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(1),((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(2));

        double vtx_x = -wcsimrootevent->GetVtx(0), vtx_y = wcsimrootevent->GetVtx(2), vtx_z = wcsimrootevent->GetVtx(1);
        // Extrapolate to get muon entrance point
        double entrance_x = vtx_x - BeamDir.x(), entrance_y = vtx_y + BeamDir.z(), entrance_z = vtx_z + BeamDir.y();
        int count = 0;
        while (!(sqrt(entrance_x*entrance_x+entrance_y*entrance_y)<max_r && fabs(entrance_z)<max_z) && count++<nSteps_max)
        {
            entrance_x += -BeamDir.x(); entrance_y += BeamDir.z(); entrance_z += BeamDir.y();
        }

        if (count>=nSteps_max) continue;

        hist_vertices->Fill(entrance_x,entrance_y,entrance_z);
        hist_vertices_xy->Fill(entrance_x,entrance_y);
        hist_vertices_yz->Fill(entrance_y,entrance_z);
        hist_vertices_zx->Fill(entrance_z,entrance_x);
    }

    TCanvas* c1 = new TCanvas();

    hist_vertices->Draw("box");
    c1->SaveAs(Form("%svertices.pdf",prefix.c_str()));

    hist_vertices_xy->Draw("colz");
    c1->SaveAs(Form("%sverticesXY.pdf",prefix.c_str()));

    hist_vertices_yz->Draw("colz");
    c1->SaveAs(Form("%sverticesYZ.pdf",prefix.c_str()));

    hist_vertices_zx->Draw("colz");
    c1->SaveAs(Form("%sverticesZX.pdf",prefix.c_str()));

    f->Close();
    t->Reset();
}