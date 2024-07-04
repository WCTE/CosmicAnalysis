# WCTE_Cosmic
Analysis code for WCTE cosmic muons. Intended to be used with WCTE software container.

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

## read_main_track.C
(Jimmy) Produces event display for a single event.
```
root -l -b -q read_main_track.C\(iEvent,\"filename\"\)
```


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