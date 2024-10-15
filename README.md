# nmr_timeseries_model
Code for "Using Machine Learning to Improve the Hard Modeling of NMR Time Series" paper

Version requirements: Matlab 2021b or later, Python 3.8.10 or later

Getting started: Run main.m to get a quick glance at the program. If you want to use machine learning, define the correct path of the virtual environment which knows torch and numpy

Loading data sets: Requires a phase and baseline-corrected data set as a .mat file. The file should contain 3 different variables.
 -- The data, called 'd' or 'D'
 -- The frequency axis, called 'x' or 'X'
 -- The time axis, called 't' or 'T'.
