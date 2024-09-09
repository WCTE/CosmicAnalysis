#define WaterCherenkovAng 0.719887812
#define WaterRefractiveIndex 1.33

//return true when the main track has entrance or exit radius larger than the cutoff, otherwise return false
//EntranceCut and ExitCut are the cutoff radius of the entrance and exit point
bool selection(const double Entrance[3], const double Exit[3], const double EntranceCut = 120, const double ExitCut = 145) {
    if ((Entrance[0]*Entrance[0]+Entrance[1]*Entrance[1])>(EntranceCut*EntranceCut)) {return true;}
    if ((Exit[0]*Exit[0]+Exit[0]*Exit[0])>(ExitCut*ExitCut)) {return true;}
    return false;
}

//
std::vector<std::vector<double>> TrackReconstruction(const double VerPos[3], const double VerDir[3], const double PMTPos[2014][3], const double PMTOri[2014][3]) {
    std::vector<std::vector<double>> PMTReTrack;
    double VerToPMTDir[3], PhoDir[3], VerToPMTDist, VerToPhoDist, PhoToPMTDist, theta, HitTime, CosInAngle;
    for (int i=0; i<2014; i++) {
        //VerToPMTDir[3] is a unit vector pointing from the given vertex to the PMT
        VerToPMTDir[0] = PMTPos[i][0] - VerPos[0];
        VerToPMTDir[1] = PMTPos[i][1] - VerPos[1];
        VerToPMTDir[2] = PMTPos[i][2] - VerPos[2];
        VerToPMTDist = sqrt(VerToPMTDir[0]*VerToPMTDir[0] + VerToPMTDir[1]*VerToPMTDir[1] + VerToPMTDir[2]*VerToPMTDir[2]);
        VerToPMTDir[0] = VerToPMTDir[0]/VerToPMTDist;
        VerToPMTDir[1] = VerToPMTDir[1]/VerToPMTDist;
        VerToPMTDir[2] = VerToPMTDir[2]/VerToPMTDist;

        //theta is the angle between the vertex direction vector and the VerToPMTDir[3] vector
        theta = VerToPMTDir[0]*VerDir[0] + VerToPMTDir[1]*VerDir[1] + VerToPMTDir[2]*VerDir[2];     //This is the value of cos(theta)
        //Reconstruction fail when cos(theta) > 1
        if (fabs(theta)>1) {
            PMTReTrack.push_back({-1});     //the array {-1} indicates failing of reconstruction
            continue;
        }
        theta = acos(theta);
        //Reconstruction fail when theta > emission angle of the Cherenkov photon
        if (theta>=WaterCherenkovAng) {
            PMTReTrack.push_back({-1});
            continue;
        }
        VerToPhoDist = sin(WaterCherenkovAng - theta)*VerToPMTDist/sin(TMath::Pi()-WaterCherenkovAng);
        PhoDir[0] = PMTPos[i][0] - (VerPos[0]+VerDir[0]*VerToPhoDist);
        PhoDir[1] = PMTPos[i][1] - (VerPos[1]+VerDir[1]*VerToPhoDist);
        PhoDir[2] = PMTPos[i][2] - (VerPos[2]+VerDir[2]*VerToPhoDist);
        PhoToPMTDist = sqrt(PhoDir[0]*PhoDir[0] + PhoDir[1]*PhoDir[1] + PhoDir[2]*PhoDir[2]);      //Distance form the starting position of the Cherenkov photon to the PMT
        PhoDir[0] /= PhoToPMTDist;
        PhoDir[1] /= PhoToPMTDist;
        PhoDir[2] /= PhoToPMTDist;
        CosInAngle = -PMTOri[i][0]*PhoDir[0]-PMTOri[i][1]*PhoDir[1]-PMTOri[i][2]*PhoDir[2];     //Cosine of the incident angle
        //ignored the case that a photon enters a PMT at the back
        if (CosInAngle<0) {
            PMTReTrack.push_back({-1});
            continue;
        }
        HitTime = (1e9)*(VerToPhoDist + PhoToPMTDist*WaterRefractiveIndex)/TMath::Ccgs();   //Hit time of the photon, where t=0 is the time when the main track start from the given vertex position
        PMTReTrack.push_back({HitTime, CosInAngle, PhoToPMTDist});
    }
    return PMTReTrack;
}