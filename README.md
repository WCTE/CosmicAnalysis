# WCTE_Cosmic
Analysis code for WCTE cosmic muons. Intended to be used with WCTE software container.

## EventDisplay_SingleEvent.c
Produces event display for a single event.
```
root -l -b -q EventDisplay_SingleEvent.c\(\"file_name\",evt_num\)
```
Example outputs are shown [here](fig/).

## VertexDistribution.c
Produces vertex distribution for all events 
```
root -l -b -q VertexDistribution.c\(\"file_name\"\)
```