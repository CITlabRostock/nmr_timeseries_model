%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spliner %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (c) 2024 Jan Hellwig
% University of Rostock, Department of Mathematics.

% The software can be used for academic, research and other similar noncommercial uses. The user
% acknowledges that this software is still in the development stage and that it is provided by the
% copyright holders and contributors "as is" and any express or implied warranties, fitness for a
% particular purpose are disclaimed. In no event shall the copyright owner or contributors be liable
% for any direct, indirect, incidental, special, exemplary, or consequential damages.
% The copyright holders provide no reassurances that the source code provided does not infringe any
% patent, copyright, or any other intellectual property rights of third parties. The copyright
% holders disclaim any liability to any recipient for claims brought against recipient by any third
% party for infringement of that parties intellectual property rights.

%% History

ver='2.0.2';

%% v2.0.2
% Version published to https://github.com/CITlabRostock/nmr_timeseries_model

% Complete overhaul of all comments
% Got rid of 'opt.pred' for doing nothing
% Fixed bug, where opt.isSaveFile did nothing
% Improved runtime of Python code

%% v2.0.1
% Implemented initial prediction of layers to optimize with ML. This is
% done after splitting the next layer into windows, squeezing the windows
% into 1000 data points, giving all splits to python, recieving a pydict
% with each window having 5 + 1 (or possibly more) sets of
% parameters for 1-5 peaks and for the extrapolated prediction and then 
% deciding which set of parameters fits best. Finally, we sort the found 
% peaks to the order in which they appear in the model

% Added variant: -1 with ML and without Multiplets, sad.
% Update: Code works with multiplets now, happy. Also: Bugfixes

% Added variant: -2 with ML and splitting the x-axis into independent
% subsystems

% Also added Wasserstein-Metric for models to deal with high and irregular
% quadratic errors.

%% v2.0.0
% Code completely rewritten from scratch
% Implementation of both Gauss- and Gauss-Lorentz-models
% Implementation of multiplets instead of many singular peaks
% Dataset was removed, data now contains value 'clim' for better viewing
% Runtime was improved, by creating lookup-table for binomial coefficients
% Minimized* vector operations in often called subroutine gl (*probably)
% Renewal of (my)sgfilt.m and (my)sg.m for better understanding
% Added option to disturb data via incorrect baseline / phase (testing purposes only)
% Re-added adapt for singlets, not so useless after all
% Added option to significantly reduce runtime via full / short use of preconditioning
% Added function to calculate the analytic integral values of each multiplet of each spectrum
% Bug in FindMinima.m fixed, where peaks smaller than epsilon were found
% Bug in GL/G/L-model.interpPara.m fixed, where centers were extrapolated out of bounds, leading to nonsensical results
% Bug in findMinima3.m fixed, where too close peaks would still appear
% Added functions to smooth data sets through time (smoothTimeXYZ in data.m)
% Added function to create residual data set (SolutionSubtraction)
% Added option to optimize (timewise) rim areas again, improves stability there

% Added variants:
% - #0:   No split optimization
% - #1-4: Added function that guesses center values quickly and added function that plots these guesses.
% - #1-4: Optimize centers first w.r.t. minimal distance to guessed centers
% - #1:   Optimize other parameters afterwards
% - #2:   Find subsytems = areas in freq. domain where multiplets are contained, optimize subsystems afterwards with centers
% - #3:   Optimize each multiplet in a moving window separately
% - #4:   Optimize each timelayer of Tsidx similar to the first; since peak positions are known, ordering is trivial

%% v1.2.4
% Bugfixes: lotsa
% Removed unnecessary bandwidths as options, adaptive window solution
% Small improvements in GUI, durability, stability, runtime
% After optimizing, new function in Spliner.plot_modlays(*,*,*) plots model in selected timelayer

%% v1.2.3
% Bugfix : Results was empty, now saves works

%% v1.2.2
% Comment changes

%% v1.2.1
% Initial layer parameters are now obtained correctly from the selected spectrum in the GUI

%% v1.2.0
% Added functions to thin out freq. / time domain
% -> thinoutX, thinoutT

%% v1.1.0
% interactive spline grid selection
% -> added SplineGridGUI(...)

%% v1.0.0
% Transition to object-oriented programming based on Spliner2
% Stepwise optimization of splines
% Prediction for intial parameters
% Windowed optimization of the last few time layer

%% Initialize paths
disp(['Spliner - Version ' ver])
warning on                              % May be switched off

addpath('class')
addpath('models')
addpath('bin')
addpath('alt_code')
addpath('data')

%% Initialize data
% Data sets are usually located in  /data.  The data-class object is 
% called 'data1'.

% For the analysis of new data sets, store data in a .mat-file with the
% data matrix 'D', the frequency axis 'X' and the time axis 'T'. I
% recommend creating a new 'load_new_dataset'-file. There, you can thin out
% both axes and cut the data to the desied size.

load_modell_kreuzkrumm5             
% load2_test_sim_2p                     
% load_ez_multiplet
% load_med_multiplet
% load_hard_multiplet
% load_grosserTest1                     % Oh no, German! :)

% data1 = besmuddle_data(1, data1);     % Possibility to add baseline and incorrect phasing to the data set. 
% data1 = smuddlephase_data(3, data1);

%% Initialize model

% mo = Lmodel(data1.T);     % Pure Lorentzmodel for simple data sets.
% mo = Gmodel(data1.T);     % You might never need a pure Gaussmodel.
mo = GLmodel(data1.T);      % Generalization of the pure Lorentzmodel, recommended.

%% Initialize data analysis object
PyPath = "/home/jan/PycharmProjects/grosserTest_Netz/venv/bin/python3"; % The path at which python3 is located inside a venv folder that knows torch and numpy.
sp = Spliner(data1,mo,PyPath);

%% Set algorithm parameters
% Default options can be found in class/Spliner.m

opt.bwdth               = 11;                 % Should equal the number of frequency points that encompasses the sharpest peak, represents window widths in freq. domain
opt.info                = 2;                  % Debug options, the bigger the number the more information is displayed, info == 2.1 / 2.2 for niche plots
opt.precon              = 1;                  % Case 1: Proper Preconditioning of cg, better when less steps desired, Case 2: No Precond. of cg, quicker steps
opt.edgeopt             = 1;                  % Boolean for optimizing first and last few layers again afterwards
opt.variant             = -1;                 % 0 for vanilla, < 0 for ML, all variants > 0 split optimization of centers and other parameters
% 1  Does nothing else
% 2  Also checks for non-overlapping areas to split optimization intervals of ppm-dimension
% 3  Optimizes every peak seperately (quick'n'dirty, not recommended)
% 4  Optimizes every layer from Tsidx similar to the first layer (quicker'n'dirtier, even less recommended)
% -1 For vanilla with initial optimizing solutions given by Neural Network
% -2 For NN initial solution and splitting of ppm-axis
opt.usemult             = 1;                  % Boolean if multiplets should be used or not. If not: Everything is modeled by single peaks

% Parameters for determining peaks and multiplets

opt.L1.d                = 3;                  % Degree of Savitzky-Golay polynomial, 0<=d<=2*M+1 is necessary, d>=3 is recommended
opt.L1.M                = opt.bwdth;          % Filter length of S-G-filter. Larger M => more noise (and maybe peaks) will be smoothed down
opt.L1.epsi             = 0.01;               % Threshold for peak height. Peaks with their maximum < epsi are ignored.
opt.L1.bw               = opt.bwdth;          % Window for searching local minimum of 2nd derivatives. There won't be peaks closer than bwdth.
opt.L1.delta            = -2e-4;              % Threshold for negativity of 2nd derivatives. Must be < 0.
opt.L1.multtol          = 5;                  % Threshold for how far the peaks of one multiplet may differ from (true) dist (number of indices)
opt.L1.verttol          = 0.05;               % Tolerance for duplet heights
opt.L1.horitol          = 0.05;               % Tolerance for multiplet distance in ppm-shift, max dist parameter
opt.L1.splitdist        = 0.2;                % Minimum distance for probable peaks to split optimization window

% Parameters for optimizing first layer / spectrum.

opt.L1opt.adapt         = 4*opt.bwdth;        % Distance of initial individual optimization window for singlets and smoothing window for center optimization (Comeback!)
opt.L1opt.amp_min       = 0;                  % Minimal amplitude of first derivative wobble for proposal of a peak (should be rather generous or 0)
opt.L1opt.amtcoeff      = 1;                  % Factor by which the number of proposed peaks is increased from the actual one, to guarantee no actual peaks are omitted
opt.L1opt.maxhw         = 1;                  % Maximum allowed halfwidth
opt.L1opt.maxcenter     = 50;                 % Maximum allowed index shift for center values to be shifted in first layer
opt.L1opt.maxdist       = opt.L1.horitol;     % Maximum allowed distance for peaks of one multiplet to be apart
opt.L1opt.maxlamb       = 0.25;               % Maximum allowed lambda value. Should not exceed 1.
opt.L1opt.minlamb       = -0.05;              % Minimum allowed lambda value. <0 allows for more general combinations. Should not exceed -.5.
opt.L1opt.optdist       = 1;                  % Boolean for whether or not to optimize dist
opt.L1opt.sshift        = 0.025;              % Optimization window width in variant 3 is equal to 2*sshift + max(centers)-min(centers)

% Parameters for large optimization routine

opt.pred_layer_count    = 2;                  % Number of previous layers that are optimized simultaneously with next layer. Should be >= 2.

opt.scheme              = 0;                  % Parameters / options given to LSQnonlin
opt.TolX                = 1e-6;               % Threshhold for convergence in parameter space
opt.TolFun              = 1e-4;               % Threshhold for convergence in function space
opt.MaxFunEvals         = 25000;              % Maximum number of evaluations of error function
opt.MaxIter             = 200;                % Maximum number of optimization steps per layer in the main routine

opt.Init.TolX           = 1e-10;              % Initial layer should be optimized with stricter criteria.
opt.Init.TolFun         = 1e-10;
opt.Init.MaxFunEvals    = 250000;
opt.Init.MaxIter        = 200;

opt.isSaveFile          = 0;                  % Bool deciding whether to save the results

if opt.variant < 0
    opt.ML.splitppm = 0.2;                    % Minimum distance between proposed peaks to split small ML-window, 2*splitppm for large splits
    opt.ML.peakadv = 0.7;                     % When calculating the quadratic error, multiply by peakadv, for the prediction with the expected number of peaks.
end

opt.epsinits = 0;                             % Only relevant for plotting

if opt.epsinits
    disp('Do you want to save initial prediction figures? Press key to continue')
    pause
end

sp.setOption(opt);                            % Initializes options

%% Starting analysis

% sp.checkData();   % Optional, checks if all spectra have same lengths

% Setting default spline grid indices
Tsidx = sp.findSplineGrid(7); % Selects some default spectra. It is recommended to change them in the next step.

% Choosing spline indices via GUI
[Tsidx, startIdx] = SplineGridGUI_bunt(Tsidx, sp.data);

% Manual selection of Tsidx and startIDX. Tsidx are indices of all selected
% spectra, startIDX is the index (in Tsidx) of the first spectrum to be 
% analyzed. Tsidx should always start with 1 and end with
% length(sp.data.T). If not, cut some spectra from data1 during creation.

% Default choice for kreuzkrumm5:
% Tsidx = [1;11;21;31;41;51;61]
% startIdx = 2;


% Obtaining initial parameters from starting spectrum
[center, ~] = sp.findInitPeakInfo(Tsidx(startIdx)); 
% Vector 'center' contains indices of found peaks.

% Obtaining multiplet parameters
if opt.usemult
    params = sp.findMultiplets(center, Tsidx(startIdx));
else
    params.centeridx = center;
    params.amt = ones(size(center));
    params.distidx = -1*ones(size(center));
end
% Struct 'params' contains 'centeridx', 'amt', 'distidx' for multiplets

% Creates lookup-table for binomial coefficients (more efficient)
sp.setBC(params.amt);
% sp.BC(i,j) = ((i-1) choose (j-1))

% Fitting multiplets on first spectrum
initLayer = sp.analyzeFirstLayer(params, Tsidx(startIdx));
% Struct 'initLayer' has similar structure to 'sp.model.init'

% Initializing sp.model with initial fit.
sp.model.initializeConst(initLayer, Tsidx, sp.options, sp.data.ma, sp.data.X);

% Possibility to edit parameters. Redommendation: Don't edit.
% sp.showpara;

%% Optimizing centers if needed

if opt.variant ~= 0
    % Creates points at which peak are probable
    sp.setCenterProbs(params.amt); % Now with added 'progress bar'
    if opt.variant > 0
        sp.runCenterOpt(startIdx);
    end
    % sp.plotCP; % Plots the probable centers
end

%% Smoothing of dataset over time dimension, generally not recommended.

% These methods try and smooth out the dataset over time, while retaining
% the structure and shape of the peaks. This is not needed and in most
% cases unhelpful.

if sp.data.timenoisy && opt.variant > 0 
    [centerpx, borders] = sp.discretizeCenters(); % 'borders' contains indices of midpoints between peaks, also edges
    %     offsetpx = sp.determineOffsets(centerpx);
    tbwdth = 101;
    sbwdth = 500;
    SG_deg = 3;

    % Several options

    %     sp.data.smoothTime(tbwdth);
    %     sp.data.smoothTimeDetailed(tbwdth, offsetpx, borders); 
    %     sp.data.smoothTimeSausage(tbwdth, sbwdth, centerpx); 
    sp.data.smoothTimeSausage2(tbwdth, sbwdth, SG_deg, centerpx); 
end

%% Variants of the main optimization routine

% Allows to add constraints to the optimization. Not recommended.
% for i=1:sp.model.lenPeak
%     sp.model.addConstraint('decrease', 100, 'Center', i);
% end
% sp.printConstraintInfo

switch opt.variant
    case 0
        % Run optimization 
        sp.runOptimizationLinPredict(startIdx);
    case 1
        % Initialize the optimization of all other parameters, then run optimization with optimized centers
        sp.initElseOpt;
        sp.runOptimizationLinPredict(startIdx);
    case 2
        % Find subsystems relevant for optimization, then optimize everything in each window
        [splits, intidx] = sp.findOptbounds;
        eps_flag = 0;
        if sp.options.epsinits
            eps_flag = 1;
            sp.options.epsinits = 0;
        end
        % Initialize the optimization routines for each subsystem
        for i=1:length(splits)-1
            disp(['Bounded Optimization #' num2str(i)])

            multidx = find(intidx == i); % Indices of multiplets present in subsystem 'i'
            if ~isempty(multidx)
                if i==length(splits)-1 && eps_flag
                    sp.options.epsinits = 1;
                end
                sp.runBoundedOpt(startIdx, splits(i), splits(i+1), multidx);
            else
                disp(['Something is wrong. No multiplet in ' num2str(i) '-th interval.'])
            end

        end
    case 3
        % Optimizing every multiplet in a certain moving interval
        for i=1:length(params.amt)

            disp(['Optimization of multiplet # ' num2str(i) ':'])

            obj = sp.runSausageOpt(startIdx, i);

        end
    case 4
        % Optimizing every spectrum for itself, interpolating the parameters will only be done to look at the final model
        TsidxTrun = Tsidx(:)'; % necessary for proper 'for'-syntax
        TsidxTrun(startIdx) = [];
        iter = 1;
        for i = TsidxTrun
            if iter == startIdx
                iter = iter+1;
            end
            disp(['Optimizing layer: ' num2str(iter)])
            sp.analyzeSingLayer(i); % Optimization process
            iter = iter+1;
        end
    case -1
        % Optimize like case 0, but the parameters for the initial guess are predicted by a Neural Network
        sp.runOptimizationMLPredict(startIdx);
    case -2
        % Find subsystems relevant for optimization, then optimize everything in each window
        [splits, multidx] = sp.findOptboundsML; % Multidx is a cell now, skips an unnecessary step
        eps_flag = 0;
        if sp.options.epsinits
            eps_flag = 1;
            sp.options.epsinits = 0;
        end
        % % Initialize the optimization routines for each subsystem
        for i=1:length(splits)-1

            disp(['Bounded ML-Optimization #' num2str(i)])
            if ~isempty(multidx{i})
                if i==length(splits)-1 && eps_flag
                    sp.options.epsinits = 1;
                end
                sp.runBoundedOptML(startIdx, splits(i), splits(i+1), multidx{i}); % Now with AI, wow
            else
                disp(['No multiplet in ' num2str(i) '-th interval.'])
            end

        end

end
%% Results and plots
sp.plot; % Plots the simplified model and data set (resource intensive)

sp.printResidues; % Prints the Frobenius norms

if ~isempty(sp.spath) % Saves the results if something exists
    save([sp.spath filesep 'results_final.mat'])
end

% sp.plot_modlay(2, 30) % Can plot single spectra with their models and a set of peaks

% sp.calcIntegrals(3, 1) % Plots (*, bool) the selected multiplet (int, *) and calculates the integral's values over time

% warning off
% data2 = sp.subModel(sp.data, 'err'); % Subtracts model from data and plots residuals
% data3 = sp.subModel(sp.data, 'tru'); % subtracts model from data w/o negatives
% warning on

% A different metric of model quality
% red_mass = (1-(data2.calcsum/data1.calcsum))*100;
% red_mass2 = (1-(data3.calcsum/data1.calcsum))*100;
% disp(['The percentage of mass reduced amounts to ' num2str(red_mass) '%'])

% Another different metric of model quality
% wmet = sp.calcwstmet();
% disp(['The Wasserstein-Metric sums to ' num2str(wmet)])
% disp(['The Wasserstein-Metric per spectrum equals ' num2str(wmet/length(data1.D))])
% disp(['The Wasserstein-Metric per spectrum on [0,1] equals ' num2str(wmet/(length(data1.D)*(data1.X{1}(end)-data1.X{1}(1))))])
