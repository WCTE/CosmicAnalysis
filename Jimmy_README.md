## read_main_track.c

Read and display rints out some basic information of the main track (muon) of a single event from a WCSim root file.

It is run by root with the following syntax: root -l -b -q read_main_track.C\(iEvent,\"filename\"\)
"filename" is the path of the WCSim source root file and iEvent is the number of the event to be read in the file.

The code first search for the main track (with track id 1), and then print out a set of basic information of the main track particle, including: track id, PDG code, PDG code of its parent particle, initial relativistic energy (in MeV), mass (in MeV/c^2), starting and stopping positions of the track in the simulation, initial direction and the information at the crossing positions of the particle with the tank boundaries. For a track passes through the detector tank, there exist at least 4 crossing with the boundary.
A 2D histogram of PMT charges is saved ("Charge_Display.pdf"), where the cylindrical array of PMTs are projected on a 2D plane. For those main track passing through the tank, the entrance (marked in red) and exit (maked in black) position will also be ploted on the charge histogram.

## multiple_events_reader.C

Read multiple WCSim root files and plot all the entrance and exit points of the main track, as well as the total charge over energy loss of the main track in the tank. 

It is run by root with the following syntax: root -l -b -q multiple_events_reader.C\(\"filename\"\)
"filename" is the path of multiple root files, by using *(Asterisk) in the input.

TChain is applied to combine all trees of same type in the given root files into one. The entrance and exit positions on the cylindrical tank are projected on as two 2D histograms, "Entrance_Position_WCSim.pdf" and "Exit_Position_WCSim.pdf" respectively. Also a 2D histogram of total PMT charges against energy loss of the main track in the tank, "Charge_vs_energy_WCSim.pdf", will be created.

## fitQun_analysis.C

Read multiple fitQun root files and plot all the extrapolated entrance and exit points form the 1-Ring fit vertices of muons.

It is run by root with the following syntax: root -l -b -q fitQun_analysis.C\(\"filename\"\)
"filename" is the path of multiple root files, by using *(Asterisk) in the input.

The vertex position and direction will be read from the sources fitQun files. The function, extrapolation(), extrapolates the vertex to the outtest boundaries of the tank to obtain the entrance and exit positions. Tracks which did not touch the tank will be eliminated by the function. Also, as the muon tracks are supposed to be approching from the top, for the targeted samples given by the simulation while writing this code, vertices with upward direction or z coordinate lower the the bottom of the tank will be filtered out. 2D histograms of the projection of entrance and exit distance will be created and saved as, "Entrance_Position_fitQun.pdf" and "Exit_Position_fitQun.pdf" respectively.

## WCSim_fitQun_preprocess.c
Read multiple WCSim root files and the corresponding fitQun root files (i.e. performing fitQun reconstruction on the WCSim results), and created a simple root file to store some necessary data of the muon tracks. It allows fast accessing to the data for later analysis purpose.

It is run by root with the following syntax: root -l -b -q WCSim_fitQun_preprocess.c\(\"fname\", \"wname\"\)
"fname" and "wname" are the paths of the fitQun and WCSim root file respectively, again using *(Asterisk) in the input. Note that the events of the WCSim and fitQun should be provided in exactly the same order, and the number of event stored in each WCSim and fitQun file is suggested to be the same. For example, the i th entry in both the j th WCSim and fitQun file given should correspond to the same event.

The code first creates a root file, "WCSim_fitQun_preprocess.root", a tree and some branches beneath it. The branches store some data related to the muon tracks, which includes fitQun vertex position, entrance point, exit point, travel direction and travel distance for both fitQun and WCSim, fitQun likelihood, total PMT charge from WCSim, as well as the angle between the fitQun and WCSim travel direction. Reconstructed event from the given fitQun file will be first read, where the function, extrapolation(), will be called. If the track is not filtered out, necessary data will be calculated and saved as an entry of the corresponding branches. Then the WCSim data of that event will be read and calculated. The tree will be finally saved in the created root file.

## s_selection.c
Read the root file created by "WCSim_fitQun_preprocess.c" and comparing the result of WCSim simulation and the corresponding fitQun reconstruction.

It is run by root with the following syntax: root -l -b -q s_selection.c\(\"pname\"\)
"pname" is the path of the root file created by "WCSim_fitQun_preprocess.c", by default pname="./WCSim_fitQun_preprocess.root".

The code first read the given file event by event. In each event the function selector(), will be called, which the main track will be selected according to the criteria (cut off) specified in the function. Multiple of plots will be created for those events being selected, including
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
