# WCTE_Cosmic
Analysis code for WCTE cosmic muons. Intended to be used with WCTE software container.

## EventDisplay_SingleEvent.c
Produces event display for a single event.
```
root -l -b -q EventDisplay_SingleEvent.c\(\"file_name\",evt_num\)
```

## VertexDistribution.c
Produces vertex distribution for all events 
```
root -l -b -q VertexDistribution.c\(\"file_name\"\)
```


## read_main_track.C
(Jimmy) Produces event display for a single event.
```
root -l -b -q read_main_track.C\(iEvent,\"filename\"\)
```


## find_max_r_and_z.C
(Jimmy) Read the geometry tree of any event and find the max_z and max_r.
```
root -l -b -q read_main_track.C\(\"filename\"\)
```