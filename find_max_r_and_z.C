R__LOAD_LIBRARY(/opt/WCSim/build/install/lib/libWCSimRoot.so)

int find_max_r_and_z(const char* filename="/work/kmtsui/wcte/cosmic/hk_flux/wcsim_mu-_0000.root")
{
    TFile* file = new TFile(filename, "read");
    TTree* geotree = (TTree*)file->Get("wcsimGeoT");
    WCSimRootGeom* geo = new WCSimRootGeom();
    geotree->SetBranchAddress("wcsimrootgeom", &geo);
    geotree->GetEntry(0);

    //find the maximum radius and z coordinate
    double max_r = 0;
    double max_z = 0;
    for (int i=0; i<geo->GetWCNumPMT(); i++) {
        WCSimRootPMT pmt;
        pmt = geo->GetPMT(i);

        double x = -pmt.GetPosition(0);
        double y = pmt.GetPosition(2);
        double r = sqrt(x*x+y*y);
        if (r>max_r) max_r = r;

        double fabs_z = fabs(pmt.GetPosition(1));
        if (fabs_z>max_z) max_z = fabs_z;
    }

    std::cout << "Maximum radius: " << max_r << "\n" <<
        "Maximum z: " << max_z << "\n";

    return 0;
}