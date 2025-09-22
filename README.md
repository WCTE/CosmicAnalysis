# ComsicAnalysis
ROOT-based analysis to study cosmic muons.

## MC
A through-going cosmic muon sample is stored at `/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/`. The production and analysis environment can be accessed by
```
/eos/experiment/wcte/MC_Production/ContainerImage/run_container.sh /eos/experiment/wcte/MC_Production/ContainerImage/softwarecontainer_workshop.sif
```
fiTQun reconstruction is run on this sample, and the output is stored in `/eos/experiment/wcte/MC_Production/v1.4.1/cosmics/fiTQun/`.

### AnalyzeCosmicsMC.c
To select muons that enter through the top cap and exit through the bottom cap, create a set of plots in `fig/` to study cut values on number of PMT hits and charge ratios in different sections (top cap, barrel, bottom cap).
```
root -l -b -q AnalyzeCosmicsMC.c
```
Then we derive a set of selection cuts:
- number of PMT hits > 1000
- Top/Total charge < 0.07
- 0.38 < Barrel/Total charge < 0.6
- 0.38 < Bottom/Total charge < 0.6

To further improve selection purity, the fiTQun reconstructed vertex and direction are extrapolated to find the muon entrance and exit positions at the tank surface. Then we only keep muons with:
- reconstructed entrance point on the top cap
- reconstructed exit point on the top cap

After the selection, 26% of all the muons remain. The selection efficiency and purity of top-down events are 85% and 93% respectively. Meanwhile over 99% of the selected events are with truth `cos(zenith_angle)>0.7`.

### convert_to_h5.py
Convert wcsim and fiTQun root files into a single h5 for ML studies. The h5 file contains the following arrays:
- `vertex`: truth vertex (nevents,3)=shape, outside the tank
- `direction`: truth direction (nevents,3)
- `entrance_pos`: truth entrance point (nevents,3)
- `exit_pos`: truth exit point (nevents,3)
- `is_top_down`: truth top-down flag (nevents)
- `momentum`: truth momentum at vertex (nevents)
- `pmtQ`: total charge per PMT (nevents,npmts)
- `pmtT`: first hit time per PMT (nevents,npmts)
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