# neuropix2023
This repository contains MATLAB codes used in XXTX1's dissertation project for NEUR0021.
We give credit to bombcell (https://github.com/Julie-Fabre/bombcell.git), braindraw (https://github.com/Julie-Fabre/braindraw.git), and brainreg (https://github.com/brainglobe/brainreg.git). Some files and lines were taken from these repositories.

This repository is purely to serve as reference for the dissertation.
- preprocessing: code used to process Neuropixels data and histology scans to basal ganglia unit spike activity.
- data_analysis: code used to produce figures, tables, and data referenced in the dissertation.
- confirmatory: requires get_experiment and get_meanfr first to respectively retrieve unit data and mean firing rates. For statistical analysis.

MATLAB R2021b was used. Add the following packages to path first to run code in this repository:
- bombcell: https://github.com/Julie-Fabre/bombcell.git
- npy-matlab: https://github.com/kwikteam/npy-matlab.git
- braindraw: https://github.com/Julie-Fabre/braindraw.git
- distributionPlot: https://uk.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m
