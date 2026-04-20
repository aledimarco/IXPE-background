# IXPE-background
Software to reject background in IXPE data following the prescription reported in the paper A. Di Marco et al., AJ 165, 143 (2023) https://iopscience.iop.org/article/10.3847/1538-3881/acba0f

After DU2 anomaly occurred on April 2025, a new version of the software has been released with improved cuts which description is reported in A. Di Marco, "The hitchhiker's guide to the IXPE data analysis", arXiv:2604.03366 (https://arxiv.org/abs/2604.03366)

## Input parameters:
- path_lv2, Input file lv2
  
- path_lv1, list of input files lv1 corresponding to the lv2 file
  
- --output, parameter values can be: rej, default value producing lv2 file including only source events, bkg producing lv2 files including only events identified as background, tag producing a new column in the lv2 file containing for each event value 1 for events identified as source and 0 for the ones identified as background

## Example:
python filter_background.py ./02008901/event_l2/ixpe02008901_det3_evt2_v01.fits ./02008901/event_l1/2008901_det3_evt1_v01_Observation.fits --output rej
