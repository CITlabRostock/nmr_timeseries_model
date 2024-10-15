%% Test data set for Training AI, (UNUSED) (3 peaks, 3 transitions)
% Have to be defined:
file='data/modell_großerTest1/grosserTest1rL.mat';

% Can be defined:
name='3p_modelDataforAI';
clim = [0,1];
comment='3 peaks, 3 transitions';

% Constructor
data1=data(file,name,clim,comment);

% If version L = large is chosen, thin out!
data1.thinoutT(1);
data1.thinoutX(10);

