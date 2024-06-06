void ReadFQ(const char * fname)
{
    TFile* f = new TFile(fname);
    TTree* t = (TTree*)f->Get("fiTQun");

    int fqnse;
    float fq1rpos[10][7][3];
    t->SetBranchAddress("fqnse",&fqnse);
    t->SetBranchAddress("fq1rpos",fq1rpos);

    for (int i=0;i<t->GetEntries();i++)
    {
        t->GetEntry(i);
        std::cout<<"fqnse = "<<fqnse<<" fq1rpos = "<<fq1rpos[0][2][0]<<" "<<fq1rpos[0][2][1]<<" "<<fq1rpos[0][2][2]<<"\n";
    }
}