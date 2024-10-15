%% Medium model data set, three crossing multiplets
% Have to be defined:
file='data/modell_multiplets_med/modell_multiplets_med.mat';

% Can be defined:
name='medMult';
clim = [0,1];
comment='3 multiplets, 4 crossings';

% Constructor
data1=data(file,name,clim,comment);

% data1.restrictX(0,4);