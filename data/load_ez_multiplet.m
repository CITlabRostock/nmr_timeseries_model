%% Simple model data set, three crossing multiplets
% Have to be defined:
file='data/modell_multiplets_ez/modell_multiplets_ez.mat';

% Can be defined:
name='ezMult';
clim = [0,1];
comment='3 multiplets, 1 crossing';

% Constructor
data1=data(file,name,clim,comment);

% data1.restrictX(1,5);