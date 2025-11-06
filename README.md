# ComsicAnalysis
ROOT-based analysis to study cosmic muons.

## MC
A through-going cosmic muon sample is stored at `/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/`. The production and analysis environment can be accessed by
```
/eos/experiment/wcte/MC_Production/ContainerImage/run_container.sh /eos/experiment/wcte/MC_Production/ContainerImage/softwarecontainer_workshop.sif
```
fiTQun reconstruction is run on this sample, and the output is stored in `/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/fiTQun/`.

### mPMT efficiency difference
To mimic the effect of reduced PMT efficiency observed in [through going beam muons](https://wcte.hyperk.ca/wg/simulation-and-analysis/meetings/2025/20250925/meeting/studying-mpmt-efficiencies-using-crossing-muon-sample/crossingmuonstudy.pdf), [MDT](https://github.com/hyperk/MDT/tree/angular_response_gain) is used to reprocess the MC data. An ad-hoc detection efficiency factor (relative to WCSim) of 45% is applied on WUT ex-situ mPMTs, and 90% otherwise. Waveform simulation and hit finding are also applied based on self-trigger conditions.

### Bad channel list
To create bad channel list per run, the `good_wcte_pmts` list in the metadata and the actual run data are used.
```
root -l -b -q CreateBadChannelList.c
```

### AnalyzeCosmicsMC.c
To select muons that enter through the top cap and exit through the bottom cap, create a set of plots in `fig/` to study cut values on number of PMT hits and charge ratios in different sections (top cap, barrel, bottom cap).
```
root -l -b -q AnalyzeCosmicsMC.c
```
Then we derive a set of selection cuts:
- number of PMT hits > 700
- Top/Total charge < 0.07
- 0.45 < Barrel/Total charge < 0.7
- 0.25 < Bottom/Total charge < 0.5

To further improve selection purity, the fiTQun reconstructed vertex and direction are extrapolated to find the muon entrance and exit positions at the tank surface. Then we only keep muons with:
- reconstructed entrance point on the top cap
- reconstructed exit point on the top cap

After the selection, 25% of all the muons remain. The selection efficiency and purity of top-down events are 68% and 79% respectively. Meanwhile over 95% of the selected events are with truth `cos(zenith_angle)>0.7`.

### convert_to_h5.py
Convert wcsim and fiTQun root files into a single h5 for ML studies. Note that the units of length, time and momentum are mm, ns and MeV/c respectively. And there is a coordinate transformation of $(x,y,z)\rightarrow (x,-z,y)$ applied.

The h5 file contains the following arrays:
- `vertex`: truth vertex (nevents,3)=shape, outside the tank
- `direction`: truth direction (nevents,3)
- `entrance_pos`: truth entrance point (nevents,3)
- `exit_pos`: truth exit point (nevents,3)
- `is_top_down`: truth top-down flag (nevents)
- `momentum`: truth momentum at vertex (nevents)
- `pmtQ`: total charge per PMT (nevents,npmts)
- `pmtT`: first hit time per PMT (nevents,npmts)
- `fq_pos`: fiTQun reconstructed vertex and time (nevents,4)
- `fq_entrance_pos`: fiTQun reconstructed entrance point (nevents,3)
- `fq_exit_pos`: fiTQun reconstructed exit point (nevents,3)
- `fq_direction`: fiTQun reconstructed direction (nevents,3)
- `fq_momentum`: fiTQun reconstructed momentum (nevents)
- `pass_selection`: pass selection flag (nevents)

To load the arrays,
```
import h5py
import torch

with h5py.File('wcte_cosmics_mc.h5','r') as f:
    entrance_pos = torch.as_tensor(f['entrance_pos'][:]) # shape=[nevents,3]
    exit_pos = torch.as_tensor(f['exit_pos'][:]) # [nevents,3]
    is_top_down = torch.as_tensor(f['is_top_down'][:]) # [nevents]
    pmtQ = torch.as_tensor(f['pmtQ'][:]) # [nevents,npmts]
    # etc.
```

## Data
To search for cosmic muons in data, a 50 ns moving time window is used with a number of hit threshold (N50) equal to 1000. Then the same charge ratio cuts (top/barrel/bottom over total charge) are applied. After that, the candidate event PMT hit information is written in a WCSim-like format for fiTQun processing.

In run 1766, which is a background run in pure water phase, there are 442k events with each event spanning 500 $\mu s$, giving a total live time of 211s. 1403 muon candidates are found. 

After fITQun processing, the reconstructed entrance/exit point cuts are applied, and 948 muon candidates remain.