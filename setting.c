void setting(double CosInAng) {
    double PMTPos[3] = {0, 0, -133.243778333};
    double PMTDir[3] = {0, 0, 1};
    double d = 25;
    double dx = d*sqrt(1-CosInAng*CosInAng), dz = d*CosInAng;
    double CenterPos[3] = {PMTPos[0]+dx, PMTPos[1], PMTPos[2]+dz};
    double n[3] = {-dx/d, 0, -dz/d};
    double rot1[3] = {n[2], 0, -n[0]};
    std::cout << "/gps/pos/rot1 " << -rot1[0] << " " << rot1[2] << " " << rot1[1] << "\n";
    std::cout << "/gps/pos/centre " << -CenterPos[0] << " " << CenterPos[2] << " " << CenterPos[1] << " cm\n";
    std::cout << "/gps/direction " << -n[0] << " " << n[2] << " " << n[1] << "\n";
}