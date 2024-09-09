#!/bin/bash

file_path="./macros/WCTE.mac"   #the file path containing the WCTE.mac file for the simulation

for i in {1..100}
    cosine=$(echo "$i*1" | awk '{printf "%.2f", $1/100}')   #cosine of the incident angle of the photon generated by the simulation w.r.t. the targeted PMT, ranging from 0 to 1 with an interval of 0.01
    root -l -q setting.c\($cosine\) > setting.txt   #call the "setting.c" macro by root to calculate the parameters for the simulation and store in "setting.txt" temporarily

    #get the parameters form the "setting.txt" file
    rot1=$(grep "/gps/pos/rot1" "./setting.txt")
    centre=$(grep "/gps/pos/centre" "./setting.txt")
    direction=$(grep "/gps/direction" "./setting.txt")

    #set the parameter in WCTE.mac for simulation
    sed -i "/\/gps\/pos\/rot1/c\\$rot1" "$file_path"
    sed -i "/\/gps\/pos\/centre/c\\$centre" "$file_path"
    sed -i "/\/gps\/direction/c\\$direction" "$file_path"

    WCSim macros/WCTE.mac macros/tuning_parameters.mac     #perform the simulation. Charge the path of the two files if neccessary
    mv wcsim.root wcsim$cosine.root     #rename the WCSim output root file by the value of the incident angle
done

rm setting.txt      #remove the auxiliary file "setting.txt"