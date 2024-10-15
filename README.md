# nmr_timeseries_model
Code for "Using Machine Learning to Improve the Hard Modeling of NMR Time Series" paper

## Requirements
MATLAB-Code:
Matlab 2021b or later, Python 3.8.10 or later

Python-Code:
Python 3.8.10 or later, torch, numpy, random, scipy

## Modeling data of NMR time series, quick start
Getting started: Run main.m to get a quick glance at the program. If you want to use machine learning, define the correct path of the virtual environment which knows torch and numpy

Loading data sets: Requires a phase and baseline-corrected data set as a .mat file. The file should contain 3 different variables.
 - The data, called `d` or `D`
 - The frequency axis, called `x` or `X`
 - The time axis, called `t` or `T`.

 Then, create new file for loading data, which can be copied and pasted from existing files and add it to file section in main.m.

 Set hyperparameters in options section (not necessary, but incorrect definition of hyperparamters leads to incorrect results).

 Then, run main in desired variant.

 Results can be saved, if `opt.isSaveFile = 1;`

## Training CNNs
New NNs can be trained using the designated Python-Code for training neural networks. Keep in mind that each new network should predict different numbers of peaks.

Once the networks are fully trained, save them under `/python/pytorch/MATLAB_net_[n]p_SD.pt`, where `[n]` is the number of peaks.
*DON'T* store unused newtorks in this folder!
