## read_main_track.c
The code print out some basic information of the main track (muon) of a single event from a WCSim root file. It is run by root with the following syntax:
root -l -b -q read_main_track.C\(iEvent,\"filename\"\)
where "filename" is the path of the WCSim source root file and iEvent is the number of the event to be read in the file.
The code first search for the main track (with track id 1), and then print out a set of basic information of the main track particle, including: track id, PDG code, PDG code of its parent particle, initial relativistic energy (in MeV), mass (in MeV/c^2), starting and stopping positions of the track in the simulation, initial direction and the information at the crossing positions of the particle with the tank boundaries. For a track passes through the detector tank, there exist at least 4 crossing with the boundary. 
A 2D histogram of PMT charges is saved, where the cylindrical array of PMTs are projected on a 2D plane. For those main track passing through the tank, the entrance (marked in red) and exit (maked in black) position will also be ploted on the charge histogram.


## multiple_events_reader.C
(Jimmy) Read multiple WCSim root file and print out the exit and entrance position, as well as the total charge against energy plot
```
root -l -b -q multiple_events_reader.C\(\"filename\"\)
```


## fitQun_analysis.C
(Jimmy) Read multiple fitQun files and print out the extrapolated exit and entrance position form the reconstructed vertex
```
root -l -b -q fitQun_analysis.C\(\"filename\"\)
```


## WCSim_fitQun.c
(jimmy) Read multiple WCSim and fitQun files and print out their diffenerce
```
root -l -b -q WCSim_fitQun.c\(\"fname\", \"wname\"\)
```

## WCSim_fitQun_preprocess.c
(jimmy) Read the required data form the WCSim and fitQun root files and store the data in a simpler root file, default file name: WCSim_fitQun_preprocess.root
```
root -l -b -q WCSim_fitQun_preprocess.c\(\"fname\", \"wname\"\)
```

## s_selection.c
(jimmy) Read the root file WCSim_fitQun_preprocess.root obtain by WCSim_fitQun_preprocess.c and comparing the result of WCSim simulation and the corresponding fitQun reconstruction. By default, pname = "./WCSim_fitQun_preprocess.root"
```
root -l -b -q s_selection.c\(\"pname\"\)
```