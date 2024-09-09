void setting(double CosInAng) {
    double PMTPos[3] = {0, 0, -133.243778333};    //Position of the centre of the targeted PMT
    double PMTOri[3] = {0, 0, 1};   //Orientation of the targeted PMT
    double d = 25;   //Distance of the centre of the circular photon source to the centre of the targeted PMT
    double dx = d*sqrt(1-CosInAng*CosInAng);    //x component of the distance d
    double dz = d*CosInAng;    //z component of the distance d
    double CenterPos[3] = {PMTPos[0]+dx, PMTPos[1], PMTPos[2]+dz};      //Position of the centre of the targeted PMT
    double n[3] = {-dx/d, 0, -dz/d};    //emission direction of the circular photon source
    double rot1[3] = {n[2], 0, -n[0]};     //The rot1 vector for determine the orientation of the circular photon source. Given by the cross product of rot2[3] cross n[3], where rot2[3]={0,1,0} by default.
    std::cout << "/gps/pos/rot1 " << -rot1[0] << " " << rot1[2] << " " << rot1[1] << "\n";
    std::cout << "/gps/pos/centre " << -CenterPos[0] << " " << CenterPos[2] << " " << CenterPos[1] << " cm\n";  
    std::cout << "/gps/direction " << -n[0] << " " << n[2] << " " << n[1] << "\n";
}