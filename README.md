# WCTE_Cosmic
Analysis code for WCTE cosmic muons. Intended to be used with WCTE software container. Plots and results are summarized [here](https://docs.google.com/presentation/d/1bXUWfKTYQnTuQvL2uMAMVh6uErd6iG4mXZFp1Xbc-1E/edit?usp=sharing) and [here](https://wcte.hyperk.ca/wg/simulation-and-analysis/meetings/2024/20240830/meeting/calibration-with-cosmic-muon/calibration_with_cosmic_muon_20240830.pdf). 

## EventDisplay_SingleEvent.c
Produces event display for a single event.
```
root -l -b -q EventDisplay_SingleEvent.c\(\"file_name\",evt_num\)
```
Example outputs are shown [here](fig/EventDisplay_SingleEvent).

## VertexDistribution.c
Muon vertices are extrapolated to reach the tank to produced the vertex distribution for all events.
```
root -l -b -q VertexDistribution.c\(\"file_name\"\)
```
Example outputs are shown [here](fig/VertexDistribution).

## read_main_track.c

Read and display rints out some basic information of the main track (muon) of a single event from a WCSim root file.

It is run by root with the following syntax:
```
root -l -b -q read_main_track.C\(iEvent,\"filename\"\)
```
"filename" is the path of the WCSim source root file and iEvent is the number of the event to be read in the file.

The code first search for the main track (with track id 1), and then print out a set of basic information of the main track particle, including: track id, PDG code, PDG code of its parent particle, initial relativistic energy (in MeV), mass (in MeV/c^2), starting and stopping positions of the track in the simulation, initial direction and the information at the crossing positions of the particle with the tank boundaries. For a track passes through the detector tank, there exist at least 4 crossing with the boundary.
A 2D histogram of PMT charges is saved ("Charge_Display.pdf"), where the cylindrical array of PMTs are projected on a 2D plane. For those main track passing through the tank, the entrance (marked in red) and exit (maked in black) position will also be ploted on the charge histogram.

## multiple_events_reader.C

Read multiple WCSim root files and plot all the entrance and exit points of the main track, as well as the total charge over energy loss of the main track in the tank. 

It is run by root with the following syntax:
```
root -l -b -q multiple_events_reader.C\(\"filename\"\)
```
"filename" is the path of multiple root files, by using *(Asterisk) in the input.

TChain is applied to combine all trees of same type in the given root files into one. The entrance and exit positions on the cylindrical tank are projected on as two 2D histograms, "Entrance_Position_WCSim.pdf" and "Exit_Position_WCSim.pdf" respectively. Also a 2D histogram of total PMT charges against energy loss of the main track in the tank, "Charge_vs_energy_WCSim.pdf", will be created.

## fitQun_analysis.C

Read multiple fitQun root files and plot all the extrapolated entrance and exit points form the 1-Ring fit vertices of muons. Details of the fiTQun variables are described [here](README_fiTQun).

It is run by root with the following syntax:
```
root -l -b -q fitQun_analysis.C\(\"filename\"\)
```
"filename" is the path of multiple root files, by using *(Asterisk) in the input.

The vertex position and direction will be read from the sources fitQun files. The function, extrapolation(), extrapolates the vertex to the outtest boundaries of the tank to obtain the entrance and exit positions. Tracks which did not touch the tank will be eliminated by the function. Also, as the muon tracks are supposed to be approching from the top, for the targeted samples given by the simulation while writing this code, vertices with upward direction or z coordinate lower the the bottom of the tank will be filtered out. 2D histograms of the projection of entrance and exit distance will be created and saved as, "Entrance_Position_fitQun.pdf" and "Exit_Position_fitQun.pdf" respectively.

## WCSim_fitQun_preprocess.c

Read multiple WCSim root files and the corresponding fitQun root files (i.e. performing fitQun reconstruction on the WCSim results), and created a simple root file storing necessary data of the muon tracks. It allows fast accessing to the data for later analysis purpose.

It is run by root with the following syntax:
```
root -l -b -q WCSim_fitQun_preprocess.c\(\"fname\", \"wname\"\)
```
"fname" and "wname" are the paths of the fitQun and WCSim root file respectively, again using *(Asterisk) in the input. Note that the events of the WCSim and fitQun must be provided in exactly the same order, and the number of event stored in each WCSim and fitQun file should be the same. For example, the i th entry in both the j th WCSim and fitQun file given should correspond to the same event.

The code first creates a root file, "WCSim_fitQun_preprocess.root", a tree and some branches beneath it. The branches store some data related to the muon tracks, which includes fitQun vertex position, entrance point, exit point, travel direction and travel distance for both fitQun and WCSim, fitQun likelihood, total PMT charge from WCSim, as well as the angle between the fitQun and WCSim travel direction. Reconstructed event from the given fitQun file will be first read, where the function, extrapolation(), will be called. If the track is not filtered out, necessary data will be calculated and saved as an entry of the corresponding branches. Then the WCSim data of that event will be read and calculated. The tree will be finally saved in the created root file.

## s_selection.c

Read the root file created by "WCSim_fitQun_preprocess.c" and comparing the result of WCSim simulation and the corresponding fitQun reconstruction.

It is run by root with the following syntax:
```
root -l -b -q s_selection.c\(\"pname\"\)
```
"pname" is the path of the root file created by "WCSim_fitQun_preprocess.c", by default pname="./WCSim_fitQun_preprocess.root".

The code first read the given file event by event. In each event the function selector(), will be called, which the main track will be selected according to the criteria (cut off) specified. There are 5 options determined by the first argument of the function, which are: 1 for top to bottom, 2 for top to bottom with cutoff, 3 for top to barrel, 4 for barrel to bottom, 5 for barrel to barrel. Multiple of plots will be created for those events being selected, including
1D histogram: 
    Distance between the entrance points of WCSim and fitQun, "entrance_diff.pdf"
    Distance between the exit points of WCSim and fitQun, "exit_diff.pdf"
    Difference of travel distances in the tank of WCSim and fitQun (i.e. fitQun distance - WCSim distance), "dist_diff.pdf"
    Total PMT charges over WCSim travel distance, "Q_over_trueL.pdf"
    Total PMT charges over fitQun travel distance, "Q_over_reconstructedL.pdf"
    Combination of the above two histograms of WCSim and fitQun, "Q_over_dist.pdf"
2D histogram:
    Total PMT charges against the WCSim travel distance, "QvsTL.pdf"
    Total PMT charges against the fitQun travel distance, "QvsRL.pdf"

## m_selection.c

Read the root file created by "WCSim_fitQun_preprocess.c" multiple times with different fitQun entrance or exit cutoff and comparing the result of WCSim simulation and the corresponding fitQun reconstruction.

It is run by root with the following syntax:
```
root -l -b -q m_selection.c\(\"pname\"\)
```
"pname" is the path of the root file created by "WCSim_fitQun_preprocess.c", by default pname="./WCSim_fitQun_preprocess.root".

The file will be read, by a for loop, by a multiple times determined by the input value of start, end and interval. Noted that (end-start) must be divisible by interval. Within the outter most for loop, the file will be read event by event and the selection() function will be called. Only main track which enters from top and exit through bottom will be selected. There are also two options determined by the first argument of the function, which are: 1 for entrance radius cutoff and 2 for exit radius cutoff. Four 1D histograms will be created for those events being selected, including:
    Distance between the WCSim and fitQun entrance position, "entrance_diff_mean.pdf"
    Distance between the WCSim and fitQun exit position, "exit_diff_mean.pdf"
    Distance between the WCSim and fitQun travel distance in the tank, "dist_diff_mean.pdf"
    Percentage of event ramaining after cutoff, "percentage.pdf"

## WCSimfitQunHitPreprocess.c

Read multiple WCSim root files and the corresponding fitQun root files (i.e. performing fitQun reconstruction on the WCSim results), and created a simple root file storing necessary data of the muon tracks as well as the true and digit hits of all PMT. It allows fast accessing to the data for later analysis purpose.

It is run by root with the following syntax:
```
root -l -b -q WCSimfitQunHitPreprocess.c\(\"fname\", \"wname\"\)
```
"fname" and "wname" are the paths of the fitQun and WCSim root file respectively, again using *(Asterisk) in the input. Note that the events of the WCSim and fitQun must be provided in exactly the same order, and the number of event stored in each WCSim and fitQun file should be the same. For example, the i th entry in both the j th WCSim and fitQun file given should correspond to the same event.

The code first creates a root file, "WCSimfitQunHitPreprocess.root", which contains 4 trees beneath it. WCSim and fitQun data of the muon tracks will be stored in the EventTree. THTree and DHTree will save the information of true and digit hit of the corresponding event while the GeoTree contains the position and orientation of all PMTs. Simular to ## WCSim_fitQun_preprocess.c, the code will read the given WCSim and the corresponding fitQun root file event by event and select the main muon track to be read according to the selection() functions. Then the ture and digit hits of that event will be read and saved.

How to access to the true and digit hit:

Each entry of the EventTree contains two branched, THNum and DHNum which constain the number of true and digit hits of this event respectively. Each true/digit hit will take up one entry of the THTree and DHTree. For example, if THNum=i and DHNum=j for the event in the 0th entry of the EventTree, the first ith entries in the THTree and the first jth entries in the DHTree will be the true/digit hits of that event, and the other follows. While the EventTree can be accessed directly at any entry, the THTree and DHTree should be accessed by, first looping the EventTree from the start and get THNum and DHNum, and use them to access the true/digit of the corresponding event.
The branch DHID is also contains in each entry of the THTree which carries the position in the DHTree of the digit hit generated by this true hit. If the true hits does not generate a digit hit, DHID will be set to -1. This provides an alternative way to access to the DHTree from the THTree. Noted that not all true hit can generate a digit hit, or in other words, looping the THTree and access to the DHTree by the DHID method method won't cover the whole DHTree.

## AngAcceptance.c

Read the root file created by "WCSimfitQunHitPreprocess.c" and plot the incident angle and angular acceptance histogram of different cases for angular acceptance analysis.

It is run by root with the following syntax:
```
root -l -b -q AngAcceptance.c\(\"filename\"\)
```
"filename" is the path of the root file created by "WCSimfitQunHitPreprocess.c", by default filename="./WCSimfitQunHitPreprocess.root"

The code first reads the given file event by event, selecting the main track by the function selection(), entrance cutoff=120cm and exit cutof=145cm by default. Next, reconstruction of the Cherenkov photon track, based on both WCSim and fitQun main track, will be performed. The expected incident angle of the Cherenkov photon of both WCSim and fitQun cases will be saved as histograms, TrueInAngle and fitQunInAngle respectively. The angular acceptance can be calculated as Q*d, where Q is the charge(number of hits of a PMT) and d is the travel distance of the Cherenkov photon to the PMT. Three types of angular acceptance will be plotted as Q*d against the cosine of incident angle and saved in the root file "AngAcceptance.root". fitQunDH is the angular acceptance with Q being the digit hit of the PMT and d is reconstructed from the fitQun maintrack. For fitQunTH, Q is the number of true hits of the PMT and d is reconstructed from the fitQun maintrack. For WCSimTH, Q is the number of true hits of the PMT and d is reconstructed from the WCSim maintrack. For the WCSimTH and fitQunDH, there is a histogram of the Q*d value for each bin (in in total, ranging form 0 to 1 with interval of 0.1), with the two case plotting together for comparison. 
Files produced by the code:
    Root file storing the three types of angular acceptance histogram, "AngAcceptance.root"
    Histogram of the expected incident angle of the Cherenkov photon from fitQun, "fitQunInAngle.pdf"
    Histogram of the expected incident angle of the Cherenkov photon from WCSim, "TrueInAngle.pdf"
    Histogram of Q*d for each angular bin for WCSimTH and fitQunDH, 10 files ranging from "BinHist[0.0-0.1].pdf" to "BinHist[0.9-1.0].pdf"

## TimeError.c

Read the root file created by "WCSimfitQunHitPreprocess.c" and plot the time error, difference between the reconstructed hit time of the Cherenkov photon and the true hit time, of both WCSim and fitQun cases.

It is run by root with the following syntax:
```
root -l -b -q TimeError.c\(\"filename\"\)
```
"filename" is the path of the root file created by "WCSimfitQunHitPreprocess.c", by default filename="./WCSimfitQunHitPreprocess.root"

The code first reads the given file, selects the main track by the function selection(), and reconstructs of the Cherenkov photon track, based on both WCSim and fitQun. The time error, reconstructed hit time - true hit time, of both case will be calculated and saved into 2 histograms.
Histograms produced by the code:
    Hit time error of the WCSim case, "fitQunTimeError.pdf"
    Hit time error of the fitQun case, "TrueTimeError.pdf"

## AngPosError.c

Read the root file created by "WCSimfitQunHitPreprocess.c" and plot the error of the incident angle, as well as the distance between the reconstructed Cherenkov photon emission position and the true emission position.

It is run by root with the following syntax:
```
root -l -b -q AngPosError.c\(\"filename\"\)
```
"filename" is the path of the root file created by "WCSimfitQunHitPreprocess.c", by default filename="./WCSimfitQunHitPreprocess.root"

The code first reads the given file, selects the main track by the function selection(), and reconstructs of the Cherenkov photon track, based on both WCSim and fitQun. The angular error, cosine of reconstructed incident angular of photon - cosine of the true incident angle, and the position error, distance between the reconstructed Cherenkov photon emission position and the true emission position, of both cases will be calculated and saved into 4 histograms.
Histograms produced by the code:
    Angular error of the WCSim case, "WCSimAngError.pdf"
    Angular error of the fitQun case, "fitQunAngError.pdf"
    Position error of the WCSim case, "WCSimPosError.pdf"
    Position error of the fitQun case, "fitQunPosError.pdf"

## setting.c

The code calculates the rot1 vector, the centre of the circular plane photon source and the emission direction of the photon for the simulation, which generates photons with the specified incident angle to the targeted PMT for angular acceptance analysis.

It is run by root with the following syntax:
```
root -l -b -q setting.c\(\"CosInAng\"\)
```
"CosInAng" is the cosine of the incident angle to the targeted PMT of the photons generated by the G4 General Particle Source.

## GenerateWCSim.sh
The shell script generates multiple simulations of photons directing to the targeted PMT, with cosine of incident angle ranging form 0 to 1 with an interval of 0.01. The output will be saved as root files named behind the consine value of the incident angle, ranging from "wcsim0.00.root" to "wcsim1.00.root". The file "setting.txt" will appear while executing the script. It will be deleted automatically when all simulations are finished. Make sure the file "setting.c" is placed in the same directory as the script. Also, a directory named "macros" should also be placed under the same working directory and the two files with path "./macros/WCTE.mac" and "./macros/tuning_parameters.mac" should be readied to start the simulations. Charge the path of the files if neccessary.

It is executed by the following syntax:
```
chmod +x GenerateWCSim.sh
./GenerateWCSim.sh

```

## ReadSim.c
The code reads all WCSim root files of the simulation of photons by produce by "GenerateWCSim.sh". Then plot the angular acceptance histogram, saved as "SimTrueAngAcceptance.pdf", with a bin width of the cosine of incident angle = 0.01. The averaged version of the histogram, with bin width 0.1, will be saved in the root file "SimAngAcceptance.root", for later comparison with the result from the cosmic muon simulation.

It is run by root with the following syntax:
```
root -l -b -q ReadSim.c
```

## Combine_plots.c
The code read 2 root files, "AngAcceptance.root" and "SimAngAcceptance.root", which contains the angular acceptance histogram of the cosmic muon simulation and direction simulation respectively. The "AngAcceptance.root" file contains "WCSimTH", "fitQunDH" and "fitQunTH", 3 histograms while "SimAngAcceptance.root" contains "AvgAngAcceptance". The code combines multiple histograms for comparison. 

It is run by root with the following syntax:
```
root -l -b -q Combine_plots.c
```

## WCTE.mac
This is the file contains the parameters for the photon simulation. 20000 photons with energy 2.5eV will be generated at each incident angle. Check the file for more information for the settings of the simulation. 