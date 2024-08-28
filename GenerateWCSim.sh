#!/bin/bash

file_path="./macros/WCTE.mac"

for i in {15..100}
do
    cosine=$(echo "$i*1" | awk '{printf "%.2f", $1/100}')
    root -l -q setting.c\($cosine\) > setting.txt

    rot1=$(grep "/gps/pos/rot1" "./setting.txt")
    centre=$(grep "/gps/pos/centre" "./setting.txt")
    direction=$(grep "/gps/direction" "./setting.txt")

    sed -i "/\/gps\/pos\/rot1/c\\$rot1" "$file_path"
    sed -i "/\/gps\/pos\/centre/c\\$centre" "$file_path"
    sed -i "/\/gps\/direction/c\\$direction" "$file_path"

    WCSim macros/WCTE.mac macros/tuning_parameters.mac
    mv wcsim.root wcsim$cosine.root
done

rm setting.txt