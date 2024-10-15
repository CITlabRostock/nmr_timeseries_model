%% Model data set (2 peaks, 0 transitions)
% Have to be defined
file='data/test_sim_2p/testdata_2peaks_krumm.mat';

% Can be defined
name='2p_modelData';
clim = [0,1];
comment='2 peaks, 0 transitions';

% Constructor
data1=data(file,name,clim,comment);

%data1.restrictT(0, 3);
