# Continuous isolated noise sources induce repeating waves in the coda of ambient noise correlations

[![DOI](https://zenodo.org/badge/601300942.svg)](https://zenodo.org/badge/latestdoi/601300942) ![python 3.11 only](https://img.shields.io/badge/python-3.11-blue)

![Beamforming results for master station IV.BRMO](./figures/Fig1_IV.BRMO.png)

This repository contains all data products, metadata, and code necessary to reproduce all figures of the manuscript "Continuous isolated noise sources induce repeating waves in the coda of ambient noise correlations" by Schippkus et al. (2023), in review.

`\manuscript` contains the revised manuscript pre-print pdf. Also available (and citable!) on [EarthArXiv](https://doi.org/10.31223/X52M20).

`\notebooks` contains three notebooks: `fig3_repeating_impulsive_source.ipynb` reproduces Figure 3, `fig4_sketch.ipynb` reproduces Figure 4, `fig6_secondary_mic roseism_stf.ipynb` reproduces Figure 6, `fig10_processing.ipynb` reproduces Figure 10, and `figs.ipynb` reproduces all other figures. Please read the instructions in the first cell of `figs.ipynb` carefully. `settings.toml` describes the parameters used for each figure. `schippkus_2023_lib.py` contains much of the logic for computing waveforms, cross-correlating them, and beamforming the cross correlations.

`\correlations` contains all cross-correlation functions our measurements are based on. These are computed from 2 years of continuous data in 2019 & 2020 between master stations `IV.BRMO` & `PL.OJC` and the Gräfenberg array `GR.GR*`.

`\figures` contains all figures as produced by the notebooks provided and used in the manuscript.

## Requirements

To run these notebooks, the following is required

* Python >= 3.11
* Scientific Python stack (scipy, matplotlib, numpy)
* obspy
* cartopy
* tqdm
* pygc (for easy great-circle computations)
* notebook

A functioning installation can be achieved, e.g., via conda by

```bash
>> conda create -n schippkus_et_al_2023 python=3.11
>> conda activate schippkus_et_al_2023
>> conda install -c conda-forge obspy cartopy tqdm pygc notebook
```
