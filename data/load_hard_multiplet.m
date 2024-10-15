%% Difficult model data set, seven crossing multiplets
% Have to be defined:
file='data/modell_multiplets_hard/modell_multiplets_hard.mat';

% Can be defined
name='hardMult';
clim = [0,1];
comment='7 multiplets, 4 crossings';

% Constructor
data1=data(file,name,clim,comment);