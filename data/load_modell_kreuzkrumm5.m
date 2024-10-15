%% Model data set (3 peaks, 2 transitions)
% Has to be defined:
file='data/modell_kreuzkrumm5/modelldata_kreuz_krumm5.mat';

% Can be defined:
name='3p_modelData';
clim = [0,1]; % Color limits, may improve visibility in the GUI.
comment='3 peaks, 2 transitions';

% Constructor
data1=data(file,name,clim,comment);

% Decreasing data size

% data1.restrictT(0, 3);
% data1.restrictX(3,9);
% data1.ThinourT(80); % Thins out time dimension to p%
% data1.ThinoutX(60); % Thins out frequency axis to p%

