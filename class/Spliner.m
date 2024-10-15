classdef Spliner < handle
    % Big class containing all relevant parameters, models and functions.
    properties
        data            % Data object
        model           % Model object
        options         % Struct containing settings for algorithm
        mpath           % Path to main file
        spath           % Path for results
        ppath           % Path of python3 in correct venv

        comptime = -1;  % Computation time
        BC              % Matrix of binomial coefficients, for saving time
    end

    methods
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Methods for all versions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Constructor
        function obj = Spliner(data, mo, pypa)
            if nargin <=2
                pypa = '';
            end
            obj.data = data;
            obj.model = mo;
            obj.BC = 0;

            % Default options, more detailed explanation in main.m. Will be
            % overwritten with the parameters defined in main.m.

            opt.bwdth               = 11;
            opt.info                = 2;
            opt.precon              = 1;
            opt.edgeopt             = 1;
            opt.variant             = -1;
            opt.usemult             = 1;

            opt.L1.d                = 3;
            opt.L1.M                = opt.bwdth;
            opt.L1.epsi             = 0.01;
            opt.L1.bw               = opt.bwdth;
            opt.L1.delta            = -1e-4;
            opt.L1.multtol          = floor(opt.bwdth/2);
            opt.L1.verttol          = 0.05;
            opt.L1.horitol          = 0.1;
            opt.L1.splitdist        = 0.2;

            opt.L1opt.adapt         = 4*opt.bwdth;
            opt.L1opt.amp_min       = 0;
            opt.L1opt.amtcoeff      = 1.2;
            opt.L1opt.maxhw         = 0.25;
            opt.L1opt.maxcenter     = 100;
            opt.L1opt.maxdist       = opt.L1.horitol;
            opt.L1opt.maxlamb       = 0.25;
            opt.L1opt.minlamb       = -0.05;
            opt.L1opt.optdist       = 1;
            opt.L1opt.sshift        = 0.025;

            opt.pred_layer_count    = 2;

            opt.scheme              = 1;
            opt.TolX                = 1e-7;
            opt.TolFun              = 1e-5;
            opt.MaxFunEvals         = 25000;
            opt.MaxIter             = 200;

            opt.Init.TolX           = 1e-12;
            opt.Init.TolFun         = 1e-12;
            opt.Init.MaxFunEvals    = 250000;
            opt.Init.MaxIter        = 500;
            opt.isSaveFile          = 1;

            opt.ML.splitppm = 0.5;
            opt.ML.peakadv = 0.9;

            opt.epsinits            = 0;

            obj.options = opt;

            obj.mpath = cd;
            obj.spath = '';
            obj.ppath = pypa; % For example, "/home/jan/PycharmProjects/MatlabNets/venv/bin/python3"
        end

        %% Setting options
        function obj =  setOption(obj,opt)
            obj.options = opt;
            obj.updateSavePath; % Since opt.isSaveFile needs to be considered
        end

        %% Setting the BC matrix
        function obj = setBC(obj, amt)
            % Creates all binomial coefficients up to m-1 choose k recursively
            m = max(amt);
            BCM = zeros(m, m);
            first = zeros(1,m);
            first(1) = 1;
            BCM(1,:) = first;
            BCM(:,1) = ones(m,1);
            for i=2:m
                for j=2:m
                    BCM(i,j) = BCM(i-1,j)+BCM(i-1,j-1);
                end
            end
            obj.BC = BCM;
        end

        %% Outputting options
        function options = getOption(obj)
            options = obj.options;
        end

        %% Updating save path and creating name for save file
        function obj = updateSavePath(obj)
            if obj.options.isSaveFile
                d = datetime('now');
                obj.spath = [cd filesep 'results' filesep obj.data.name '_' num2str(year(d)) '_' num2str(month(d)) '_' num2str(day(d)) '_' num2str(hour(d)) '_' num2str(minute(d)) '_' num2str(round(second(d)))];
                warning off
                mkdir([ cd filesep 'results'])
                warning on
                mkdir(obj.spath);
            else
                obj.spath='';
            end
        end

        %% Checking data
        function checkData(obj)
            if obj.data.check()
                pause(1)
            else
                pause(5) % Something's wrong.
            end
        end

        %% Setting up and calculating peak (not multiplet) centers of initial spectrum
        function [center, specidx] = findInitPeakInfo(obj, specidx_inp)

            disp('### Determine initial peak information (number of peaks, centers) ###')
            try
                disp(['# Spectrum ' num2str(specidx_inp) ' was chosen: '])
            catch
            end
            d       = obj.options.L1.d;
            M       = obj.options.L1.M;
            N       = 2*M+1;
            epsi    = obj.options.L1.epsi;
            bw      = obj.options.L1.bw;
            delta   = obj.options.L1.delta;

            if nargin == 1
                disp('Please define an initial spectrum, defaulting to first spectrum.')
                specidx_inp = 1;
            end
            Xinit   = obj.data.X{specidx_inp};
            Dinit   = obj.data.D{specidx_inp};

            center = findStartSpec4(Xinit, Dinit, bw, epsi, delta, d, N, obj.options.info); % Last argument for plotting

            specidx = specidx_inp;
            lenPeak = length(center);
            disp(['# Number of peaks = ' num2str(lenPeak)])
        end

        %% Finding out information about multiplets
        function params = findMultiplets(obj, center, specidx)
            % Center contains ascending values of peak centers
            Dinit = obj.data.D{specidx};
            Xinit = obj.data.X{specidx};
            bw = obj.options.bwdth;
            tol = obj.options.L1.multtol;
            vtol = obj.options.L1.verttol;
            htol = obj.options.L1.horitol;
            params.centeridx = [];
            params.amt = [];
            params.distidx = [];
            multiplets = zeros(length(center)); % Only used for debugging purposes
            unassigned = center; % Remove center as soon as it's been known to be part of a multiplet
            unchecked = center; % Remove center as soon as it's been checked

            % Gathering noise parameters
            addcenter = [0 center length(Xinit)];
            distances = diff(addcenter);
            noi = find(distances == max(distances),1);
            middle = ceil((addcenter(noi+1)+addcenter(noi))/2);
            expvalue = mean(Dinit(middle-bw:middle+bw)); % Minimum-variance unbiased estimator for expectation value
            variance = 1/(2*bw)*(sum((Dinit(middle-bw:middle+bw)-expvalue).^2)); % M.-v.u.e. for variance
            stddevia = sqrt(variance);


            % Is peak i part of a (3+)-let?
            % Needs adjustment parameter for distances to still be part of same multiplet
            % also checks if multiplets have ascending, then descending peak heights,
            % otherwise a false positive was found (meaning: distances fit coincidentally)
            actcenter = center; % Only has unassigned peaks in it
            while ~isempty(unchecked)
                peak = unchecked(1);
                polyplet = peak;
                ind = find(peak == center);
                if ind >= length(center)-1 % Less then 3 left to end of center
                    break;
                end
                for j = ind+1:length(center)-1
                    if abs(Xinit(center(j))-Xinit(center(ind))) > htol % Otherwise we run the risk of finding peaks which are too far apart
                        break; % Upcoming centers are even further away
                    end
                    if actcenter(j) ~= 0 % We don't want peaks that are already chosen
                        testdist = center(j)-center(ind);
                        found = j;
                        while ~isempty(found)
                            polyplet = [polyplet center(found)]; %#ok<AGROW>
                            temp = actcenter - center(found) - testdist; % Peaks we've passed have negative values in temp
                            minval = min(abs(temp(temp > -1*testdist))); % Only look for the nearest value to the right of center(found)
                            found = [];
                            if minval <= tol
                                found = find(abs(temp) == minval, 1, 'first'); % Similar to argmin(temp)
                            end
                            testdist = center(found)-polyplet(end); % If not, errors may amount to something bigger than 'tol'
                        end
                        if length(polyplet) > 2
                            isayso = 1;
                            while(isayso) % Check if these have multiplet's form, (asc. then desc.) if not, remove last peak and try again
                                isayso = 0;
                                heights = Dinit(polyplet);
                                hdiffs = diff(heights);
                                b = 1;
                                for i=1:ceil(length(polyplet)/2)-1
                                    if hdiffs(i) + 3*stddevia < 0 || hdiffs(end-i+1)-3*stddevia  > 0 % Bonus noise deviation
                                        b = 0; % Not in multiplet's form
                                        if length(polyplet) > 3 % Remove last peak and try again
                                            polyplet = polyplet(1:end-1);
                                            isayso = 1;
                                        end
                                        break;
                                    end
                                end
                            end
                            if b
                                for i=polyplet
                                    unchecked = unchecked(unchecked ~= i); % Gets rid of found peaks
                                    unassigned = unassigned(unassigned ~= i);
                                end
                                for k=1:length(center)
                                    if isempty(find(center(k) == unassigned,1))
                                        actcenter(k) = 0; % Marks that some centers have been found
                                    end
                                end
                                params.centeridx = [params.centeridx peak]; % Define parameters of found multiplets
                                params.amt = [params.amt length(polyplet)];
                                params.distidx = [params.distidx round(mean(abs(diff(polyplet))))]; % Average distance, rounded to indices
                                multiplets(length(params.amt),1:length(polyplet)) = polyplet;
                                break; % Stops inner for-loop
                            end
                        end
                        polyplet = peak;
                    end
                end
                unchecked = unchecked(unchecked~=peak);
            end


            % Is i peak of a duplet?
            unchecked = unassigned;
            while ~isempty(unchecked)
                peak = unchecked(1);
                ind = find(peak == center);
                temp = zeros(size(actcenter));
                for i=1:length(actcenter)
                    if actcenter(i) == 0 || i <= ind || abs(Xinit(actcenter(i))-Xinit(peak)) > htol % If peaks of similar height are too far apart, ignore
                        temp(i) = 10; % Data matrices are scaled below 1
                    else
                        temp(i) = abs(Dinit(peak)-Dinit(center(i)));
                    end
                end
                found = find(temp <= 3*stddevia + vtol,1); % Nearest index, duplets are found from left to right ascending (indexwise)
                if ~isempty(found)
                    unchecked = unchecked(unchecked ~= peak);
                    unchecked = unchecked(unchecked ~= center(found));
                    unassigned = unassigned(unassigned ~= peak);
                    unassigned = unassigned(unassigned ~= center(found));
                    for k=1:length(center)
                        if isempty(find(center(k) == unassigned,1))
                            actcenter(k) = 0;
                        end
                    end
                    params.centeridx = [params.centeridx peak];
                    params.amt = [params.amt 2];
                    params.distidx = [params.distidx center(found)-peak];
                    multiplets(length(params.amt),1:2) = [peak, center(found)];
                else
                    unchecked = unchecked(unchecked ~= peak);
                end
            end

            % Other peaks without a family :(
            for peak=unassigned
                params.centeridx = [params.centeridx peak];
                params.amt = [params.amt 1];
                params.distidx = [params.distidx -1];
                multiplets(length(params.amt), 1) = peak;
            end

            [sorted, order] = sort(params.centeridx); % Sort centers in ascending order
            params.centeridx  = sorted;
            params.amt = params.amt(order);           % Sorts in same order
            params.distidx = params.distidx(order);

            disp(['# Number of multiplets = ' num2str(length(params.amt))])
        end

        %% Finding initial spline layers
        function Tsidx = findSplineGrid(obj, numb)
            Tsidx = round(linspace(1,length(obj.data.T), numb));
        end

        %% Calling first-layer-analysis function of current model
        function initLayer = analyzeFirstLayer(obj, params, specidx)
            Xinit = obj.data.X{specidx}(:); %might not be necessary, will solve while constructing data
            Dinit = obj.data.D{specidx}(:);
            % no need to save in spath when it will be overwritten anyways
            initLayer = obj.model.analyzeFirstLayer(params, Xinit, Dinit, obj.options, [], obj.BC);
        end

        %% Calling show-object function of current model, allows manual change of parameters
        function obj = showpara(obj)
            currver = version('-release');
            if str2double(currver(1:4)) < 2018
                warning('This feature requires Matlab 2018 or newer')
                return
            end
            obj = obj.model.showpara;
        end

        %% Calling optimizing-whole-model function for all of the model parameters at once, UNUSED, not recommended.
        function runOptimization(obj)
            obj.comptime = obj.model.optimizeModel(obj.options, obj.data.X, obj.data.D, obj.data.ma, obj.spath, obj.BC);
        end

        %% Big optimization routine (only few time layers at a time are optimized)
        function obj = runOptimizationLinPredict(obj, TsStartIdx)
            % TsStartIdx is the Index of Tsidx to start with
            opt = obj.options;
            X = obj.data.X;
            D = obj.data.D;
            ma = obj.data.ma;
            Tsidx = obj.model.Tsidx;

            % 'SpStatus': 0 ~ unoptimized, 1 ~ upcoming, 2 ~ already optimized
            SpStatus = zeros(size(Tsidx));
            SpStatus(TsStartIdx) = 1;

            % Running initial optimization
            disp(['Run initial: Selected layer # ' num2str(TsStartIdx)])
            optIDX = TsStartIdx;
            tempMod = obj.model.getReducedModel(optIDX);
            tempMod.init = tempMod.opt;
            opttime = tempMod.optimizeModel(opt, X(Tsidx(optIDX)), D(Tsidx(optIDX)), ma(Tsidx(optIDX)), obj.spath, obj.BC); % 'X(inds)' etc. gives cell array, not its contents
            disp(['Computation time: ' num2str(opttime) 's'])
            % Transferring results back to original model
            obj.model.insertPara(tempMod.opt, optIDX);
            SpStatus(optIDX) = 2;

            % Running forward
            for j = TsStartIdx +1 : length(Tsidx)
                disp(['Run forward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = max(TsStartIdx,j-opt.pred_layer_count);
                ide = j;
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getReducedModel(optIDX);
                tempMod.init = tempMod.opt;
                xmin = X{1}(1);
                xmax = X{1}(end);
                tempind = find(SpStatus(optIDX)==1);
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax); % Interpolates Parameters

                % Plots
                if opt.info == 2.2
                    disp(['Plot initial prediction of layer # ' num2str(Tsidx(j))])
                    d = D{Tsidx(j)};
                    x = X{Tsidx(j)};
                    MLmodel = zeros(size(x));
                    for p = 1:size(tempMod.init.center, 1)
                        MLmodel = MLmodel + glmult(x, tempMod.init.amt(p, tempind), tempMod.init.center(p, tempind), tempMod.init.dist(p, tempind), tempMod.init.hw(p, tempind), tempMod.init.height(p, tempind), tempMod.init.lambda(p, tempind), obj.BC);
                    end
                    figure(7);
                    set(gcf, 'Position', get(0, 'Screensize'));
                    if opt.epsinits
                        clf
                        hold on
                        plot(x,d,'k', 'LineWidth', 2)
                        plot(x,MLmodel, 'b', 'LineWidth', 2)
                    else
                        plot(x, d, 'k', x, MLmodel, 'b');
                    end
                    set(gca, 'XDir', 'reverse')
                    set(gcf, 'renderer', 'painters');
                    title('t1')
                    xlabel('x1')
                    ylabel('y1')
                    set(gca, 'Fontsize', 24)
                    if opt.epsinits
                        ymax = max(max(d),max(MLmodel));
                        ylim([0,ymax])
                        set(gca, 'Fontsize', 48)
                        path = ['Bilder/Inits/abc123', num2str(Tsidx(j))];
                        pause(2)
                        print('-depsc', '-loose', path);
                        hold off
                    end
                    wm = wasserstein(x, MLmodel, d)*100/(x(end)-x(1)); % Normalized Wasserstein metric for one spectrum
                    frobm = norm(MLmodel - d, 'fro')*100/length(x); % Mean squared error
                    disp(['Wasserstein-Metric of: ' num2str(wm) ' %'])
                    disp(['MSE of:' num2str(frobm) ' %'])
                end

                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertPara(tempMod.opt, optIDX);
                SpStatus(j) = 2;
            end

            % Running backward
            for j = TsStartIdx -1:-1:1
                disp(['Run backward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = j;
                ide = min(length(SpStatus),j+opt.pred_layer_count);
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getReducedModel(optIDX);
                tempMod.init = tempMod.opt;
                xmin = X{1}(1);
                xmax = X{1}(end);
                tempind = find(SpStatus(optIDX)==1);
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax);

                % Plots
                if opt.info == 2.2
                    disp(['Plot initial prediction of layer # ' num2str(Tsidx(j))])
                    d = D{Tsidx(j)};
                    x = X{Tsidx(j)};
                    MLmodel = zeros(size(x));
                    for p = 1:size(tempMod.init.center, 1)
                        MLmodel = MLmodel + glmult(x, tempMod.init.amt(p, tempind), tempMod.init.center(p, tempind), tempMod.init.dist(p, tempind), tempMod.init.hw(p, tempind), tempMod.init.height(p, tempind), tempMod.init.lambda(p, tempind), obj.BC);
                    end
                    figure(7);
                    set(gcf, 'Position', get(0, 'Screensize'));
                    if opt.epsinits
                        clf
                        hold on
                        plot(x,d,'k', 'LineWidth', 2)
                        plot(x,MLmodel, 'b', 'LineWidth', 2)
                    else
                        plot(x, d, 'k', x, MLmodel, 'b');
                    end
                    set(gca, 'XDir', 'reverse')
                    set(gcf, 'renderer', 'painters');
                    title('t1')
                    xlabel('x1')
                    ylabel('y1')
                    set(gca, 'Fontsize', 24)
                    if opt.epsinits
                        ymax = max(max(d),max(MLmodel));
                        ylim([0,ymax])
                        set(gca, 'Fontsize', 48)
                        path = ['Bilder/Inits/abc124', num2str(Tsidx(j))];
                        pause(2)
                        print('-depsc', '-loose', path);
                        hold off
                    end
                    wm = wasserstein(x, MLmodel, d)*100/(x(end)-x(1)); % Normalized Wasserstein metric for one spectrum
                    frobm = norm(MLmodel - d, 'fro')*100/length(x); % MSE
                    disp(['Wasserstein-Metric of: ' num2str(wm) ' %'])
                    disp(['MSE of:' num2str(frobm) ' %'])
                end

                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertPara(tempMod.opt, optIDX);
                SpStatus(j) = 2;
            end

            if obj.options.edgeopt && opt.pred_layer_count < length(Tsidx) % Optimize first and last few layers extra
                for j=2:opt.pred_layer_count
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(1:j)])

                    optIDX = 1:j;

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getReducedModel(optIDX);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertPara(tempMod.opt, optIDX);
                end
                for j=length(Tsidx)-1:-1:length(Tsidx)-opt.pred_layer_count+1
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(j:length(Tsidx))])

                    optIDX = j:length(Tsidx);

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getReducedModel(optIDX);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertPara(tempMod.opt, optIDX);
                end
            end
        end

        %% Calling plot function of current model
        function obj = plot(obj)
            obj.model.plot(obj); % We indeed need to give 'model' the Spliner-object
        end

        %% Calling the print-constraint-info function of current model
        function printConstraintInfo(obj)
            obj.model.printConstraintInfo;
        end

        %% Plotting chosen multiplets and data at a selected timelayer, highlighting one multiplet throughout
        function obj = plot_modlay(obj, peakno, Tidx, peakind)
            if nargin == 3
                peakind = 1:obj.model.lenPeak;
            end
            peakind = peakind(:)'; % Guarantees a row vector
            Tsidx = obj.model.Tsidx;
            k = Tidx;
            try
                Xakt = obj.data.X{k};
                Dakt = obj.data.D{k};
            catch
                disp('Index of timelayer exceeds bounds.')
                return;
            end
            if length(Tidx)~= 1
                disp('Only one timelayer allowed. Why not keep seperate figures open? :p')
                return;
            end
            figure();
            col_1 = 'b';
            if obj.options.info == 2.1 % Niche plots
                hold on
                col_1 = 'k';
            end
            set(gca, 'XDir', 'reverse');
            plot(Xakt, Dakt, 'color', col_1);
            hold on
            Names = obj.model.ParaName;
            if ~isempty(find(k == Tsidx, 1)) % No need to interpolate
                for j = peakind
                    i = find(k == Tsidx); %#ok<NASGU>
                    for ind = 1:length(Names)
                        eval([Names{ind} '1 = obj.model.opt.' Names{ind} '(j,i);'])
                    end
                    switch obj.model.mname % Lambda is already set by loop above in Gauss-Lorentzmodel
                        case 'Gaussmodel'
                            lambda1 = 1;
                        case 'Lorentzmodel'
                            lambda1 = 0;
                    end
                    y = glmult(Xakt, amt1, center1, dist1, hw1, height1, lambda1, obj.BC);
                    if j == peakno
                        plot(Xakt, y, 'color', 'g');
                    else
                        plot(Xakt, y, '--',  'color', 'r');
                    end
                end
            else % Need to interpolate params with pchip
                T = obj.data.T; %#ok<NASGU>
                for j = peakind
                    for ind = 1:length(Names)
                        eval([Names{ind} '_all = pchip(T(Tsidx), obj.model.opt.' Names{ind} '(j,:), T);'])
                    end
                    switch obj.model.mname
                        case 'Gaussmodel'
                            lambda_all = ones(size(hw_all));
                        case 'Lorentzmodel'
                            lambda_all = zeros(size(hw_all));
                    end
                    y = glmult(Xakt, amt_all(Tidx), center_all(Tidx), dist_all(Tidx), hw_all(Tidx), height_all(Tidx), lambda_all(Tidx), obj.BC);
                    if j == peakno
                        plot(Xakt, y, 'color', 'g');
                    else
                        plot(Xakt, y, '--', 'color', 'r');
                    end
                end
            end
            set(gca, 'XDir', 'reverse')
            axis tight
            hold off
        end

        %% Calculating and plotting integral values of modeled multiplets, highlighting one integral throughout
        function Integrals = calcIntegrals(obj, peakno, plotcheck, peakind)
            if nargin < 4
                peakind = 1:obj.model.lenPeak;
            end
            peakind = peakind(:)'; % Guarantees a row vector
            if ~find(peakno == peakind)
                disp('Second argument must be contained in fourth argument.')
            end
            T = obj.data.T;
            Tsidx = obj.model.Tsidx; %#ok<NASGU>
            Integrals = zeros(length(T), obj.model.lenPeak);
            Names = obj.model.ParaName;
            scalefac = obj.data.scal; % For proper integrals scale back to original data
            for j = peakind
                for ind = 1:length(Names)
                    eval([Names{ind} '_all = pchip(T(Tsidx), obj.model.opt.' Names{ind} '(j,:), T);']) % Get all parameters
                end
                switch obj.model.mname
                    case 'Gaussmodel'
                        lambda_all = ones(size(hw_all));
                    case 'Lorentzmodel'
                        lambda_all = zeros(size(hw_all));
                end
                amt = amt_all(1);
                temp = height_all.*hw_all.*(((1-lambda_all).*pi)+lambda_all.*sqrt(pi/log(2))); % Analytical calculation
                Integrals(:, j) = temp./obj.BC(amt, ceil(amt/2)).*2^(amt-1).*scalefac;
            end
            % Plots
            if plotcheck
                figure();
                hold on
                suum = zeros(length(T),1);
                for j = peakind
                    suum = suum + Integrals(:,j);
                    if j == peakno
                        plot(T, Integrals(:,j), 'color', 'g', 'LineWidth',2);
                    else
                        plot(T, Integrals(:,j), 'color', 'r', 'LineWidth', 2);
                    end
                end
                plot(T, suum, 'color', 'b', 'LineWidth',2);
                xlabel('Time')
                ylabel('Peak Integrals')
            end
        end

        %% Calculating and printing residues of modeled values
        function printResidues(obj)
            para = obj.model.opt; %#ok<NASGU>
            X = obj.data.X;
            obj.model.scaleM = eye(obj.model.lenPara);

            % Reshaping and reordering of parameters
            div = obj.model.lenPara/length(obj.model.ParaName); %#ok<NASGU>
            for ind = 1:length(obj.model.ParaName)
                eval(['para_' obj.model.ParaName{ind} ' = reshape(para.' obj.model.ParaName{ind} ',[div 1]);'])
            end
            if strcmp(obj.model.mname,'Gauss-Lorentzmodel') % Compares strings
                para_2 = [para_amt; para_center; para_dist; para_hw; para_height; para_lambda];
            else
                para_2 = [para_amt; para_center; para_dist; para_hw; para_height];
            end
            % Modelling data
            M = obj.model.MOD(para_2, X, obj.BC);

            ma = obj.data.ma;
            D = obj.data.D;
            leng = 0;
            for j = 1:length(X)
                leng = leng+length(X{j}); %Just in case it's not rectangular
            end
            R_wo_scale  = zeros(leng,1);
            R_w_scale   = zeros(leng,1);
            Dref_wo     = zeros(leng,1);
            Dref_w      = zeros(leng,1);
            ir          = 0;
            for j=1:length(X)
                R_w_scale(ir+1:ir+length(X{j})) = 1/ma(j)*(D{j}-M{j});
                R_wo_scale(ir+1:ir+length(X{j})) = D{j}-M{j};
                Dref_w(ir+1:ir+length(X{j})) = 1/ma(j)*D{j};
                Dref_wo(ir+1:ir+length(X{j})) = D{j};
                ir = ir+length(X{j});
            end
            R_wo = norm(R_wo_scale, 'fro')/norm(Dref_wo,'fro');
            R_w  = norm(R_w_scale, 'fro')/norm(Dref_w,'fro');
            disp('Relative residue of unscaled matrix (w.r.t. the actual Frobenius norm):')
            disp(R_wo)
            disp('Relative residue of scaled matrix:')
            disp(R_w)
        end

        %% Subtracting found solvents
        function data2 = subModel(obj, data, mode)
            if nargin <= 1
                data = obj.data;
            end
            if nargin <= 2
                mode = 'err';
            end
            data2 = data.copy(); % Is necessary to preserve changes in classes inheriting from handle
            para = obj.model.opt; %#ok<NASGU>
            X = obj.data.X;
            obj.model.scaleM = eye(obj.model.lenPara);

            % Reshaping and reordering of parameters
            div = obj.model.lenPara/length(obj.model.ParaName); %#ok<NASGU>
            for ind = 1:length(obj.model.ParaName)
                eval(['para_' obj.model.ParaName{ind} ' = reshape(para.' obj.model.ParaName{ind} ',[div 1]);'])
            end
            if strcmp(obj.model.mname,'Gauss-Lorentzmodel') % Compares strings
                para_2 = [para_amt; para_center; para_dist; para_hw; para_height; para_lambda];
            else
                para_2 = [para_amt; para_center; para_dist; para_hw; para_height];
            end
            % Modelling data
            M = obj.model.MOD(para_2, X, obj.BC);
            for i=1:length(M)
                switch mode
                    case 'err'
                        data2.D{i} = abs(data.D{i}-M{i}); % With abs better to view?
                    case 'slv'
                        data2.D{i} = max(0, data.D{i}-M{i});
                    case 'tru'
                        data2.D{i} = data.D{i}-M{i};
                end
            end
            data2.plot();
        end

        %% Calculating the Wasserstein-metric for whole data set, might be preferrable to sum of squares
        function wmet = calcwstmet(obj)
            para = obj.model.opt; %#ok<NASGU>
            X = obj.data.X;
            obj.model.scaleM = eye(obj.model.lenPara);

            % Reshaping and reordering of parameters
            div = obj.model.lenPara/length(obj.model.ParaName); %#ok<NASGU>
            for ind = 1:length(obj.model.ParaName)
                eval(['para_' obj.model.ParaName{ind} ' = reshape(para.' obj.model.ParaName{ind} ',[div 1]);'])
            end
            if strcmp(obj.model.mname,'Gauss-Lorentzmodel') % Compares strings
                para_2 = [para_amt; para_center; para_dist; para_hw; para_height; para_lambda];
            else
                para_2 = [para_amt; para_center; para_dist; para_hw; para_height];
            end
            % Modelling data
            M = obj.model.MOD(para_2, X, obj.BC);

            % Calculating wasserstein-metric
            D = obj.data.D;
            wmet = 0;
            for i=1:length(M)
                wmet = wmet + wasserstein(X{i}, D{i}, M{i});
            end
        end
        %% Calculating the indices of the peaks and calculating borders between multiplets
        function [centerpx, borders] = discretizeCenters(obj)
            lp = obj.model.lenPeak;
            numpx = length(obj.data.X{1});
            numtl = length(obj.data.X);
            centerpx = zeros(numtl,lp);
            borders = zeros(numtl,lp+1);
            borders(:,1) = ones(numtl,1);
            borders(:,end) = numpx*ones(numtl,1);
            X = obj.data.X{1};
            amts = obj.model.opt.amt(:,1);
            for i=1:lp
                centerfull = pchip(obj.model.Ts, obj.model.opt.center(i,:), obj.model.Tfull);
                if amts(i) > 1
                    distfull = pchip(obj.model.Ts, obj.model.opt.dist(i,:), obj.model.Tfull);
                    centerfull = centerfull+((amt-1)/2)*distfull;
                end
                for j=1:numtl
                    centerpx(j,i) = find(X >= centerfull(j),1); % This does cause additional error
                    if abs(X(centerpx(j,i))-centerfull(j))>abs(X(max(centerpx(j,i)-1,1))-centerfull(j)) % Repair above error
                        centerpx(j,i) = centerpx(j,i)-1;
                    end
                end
            end
            for jj = 1:numtl
                currcpx = sort(centerpx(jj,:)); % Borders are always inbetween the i-th and (i+1)-th peak per spectrum
                borders(jj,2:lp) = floor(currcpx(1:end-1)+diff(currcpx)/2);
            end
        end

        %% Calculating offsets of each peaks sphere of influence
        function offsetpx = determineOffsets(~, cpx)
            offsetpx = zeros(size(cpx));
            numtl = size(cpx,1);
            numpks = size(cpx,2);
            cpx = sort(cpx,2); % So only the j-th-smallest peak from the left is taken into account
            for j=1:numpks
                for i=1:numtl
                    offsetpx(i,j) = cpx(i,j)-cpx(1,j);
                end
            end
        end

        %% Calculating correcting factors of noisy data
        function factors = calcFactors(obj, sbwdth, centerpx)
            numtl = length(obj.data.T);
            numpx = length(obj.data.D{1});
            tempfacs = ones(numtl,1);
            for i=1:obj.model.lenPeak
                for j=-sbwdth:sbwdth
                    offsetcol = zeros(1,numtl);
                    for t=1:numtl
                        ind = max(1,min(numpx, centerpx(t,i)+j));
                        offsetcol(t) = obj.data.D{t}(ind);
                    end
                    smoothcol = mysgfilt(1,41,offsetcol');
                    divvec = smoothcol./offsetcol';
                    divvec(divvec > 5) = 1; % For better conditioning
                    divvec(divvec < 0.2) = 1;
                    tempfacs = tempfacs + divvec;
                end
            end
            factors = tempfacs./(obj.model.lenPeak*(2*sbwdth+1));
        end

        %% Plots column (might be skewed originally)
        function plotskewcol(obj, inds)
            v = obj.data.calcskewcol(inds);
            figure();
            plot(obj.data.T,v,'b');
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variants (various) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Defining the Matrix for 'center' optimization
        function obj = setCenterProbs(obj, amt)
            disp('Creating approximates for peak positions. Completion:')
            CP_temp = obj.data.CP;
            tempidx = CP_temp;
            X = obj.data.X;
            opt = obj.options;
            nop = ceil(sum(amt)*opt.L1opt.amtcoeff); % Number of (sub-)peaks it should guess
            bw = opt.bwdth;
            d = opt.L1.d;
            amp_min = opt.L1opt.amp_min;
            warning('off') % Due to too large numbers, matrices are ill-conditioned. :(  We'll ignore that for now  :)
            prog = zeros(1,10);
            for j = 1:10
                prog(j) = ceil(length(obj.data.T)*j/10);
            end
            for i = 1:length(obj.data.T)
                [~, sg_d1, ~, ~] = mysgfilt(d,2*bw+1,obj.data.D{i});
                tempidx{i} = prob_test3(obj.data.D{i}, sg_d1, bw, nop, amp_min);
                CP_temp{i} = X{i}(tempidx{i}); % Real ppm-values of proposed peaks
                if find(prog == i)
                    disp([num2str(find(prog==i)*10) '%'])
                end
            end
            warning('on')
            if opt.info >= 3
                figure()

                hold on
                for j = 1:length(tempidx)
                    plot(tempidx{j},j*ones(size(tempidx{j})), 'rx')
                end
                set(gca, 'XDir', 'reverse');
                axis tight
            end
            obj.data.CP = CP_temp;
        end

        %% Optimization for 'center' and 'dist' values
        function obj = runCenterOpt(obj, TsStartIdx)
            % TsStartIdx is the Index of Tsidx to start with
            opt = obj.options;
            X = obj.data.X;
            CP = obj.data.CP;
            Tsidx = obj.model.Tsidx;

            % SpStatus: 0 ~ unoptimized, 1 ~ upcoming, 2 ~ already optimized
            SpStatus = zeros(size(Tsidx));
            SpStatus(TsStartIdx) = 1;

            % Running initial optimization
            disp(['Run initial: Selected layer # ' num2str(TsStartIdx)])
            optIDX = TsStartIdx;
            tempMod = obj.model.getReducedModel(optIDX);
            tempMod.init = tempMod.opt;
            opttime = tempMod.optimizeCenters(opt, X(Tsidx(optIDX)), CP(Tsidx(optIDX))); % 'X(inds)' etc. gives cell array, not its contents
            disp(['Computation time: ' num2str(opttime) 's'])
            % Transferring results back to original model
            obj.model.insertPara(tempMod.opt, optIDX); % Harmless, only uses model.opt
            SpStatus(optIDX) = 2;

            % Running forward
            for j = TsStartIdx +1 : length(Tsidx)
                disp(['Run forward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = max(TsStartIdx,j-opt.pred_layer_count);
                ide = j;
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getReducedModel(optIDX);
                tempMod.init = tempMod.opt;
                xmin = X{1}(1);
                xmax = X{1}(end);
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax);  % Also Harmless

                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeCenters(opt, X(range), CP(range));
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertPara(tempMod.opt, optIDX);
                SpStatus(j) = 2;
            end

            % Running backward
            for j = TsStartIdx -1:-1:1
                disp(['Run backward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = j;
                ide = min(length(SpStatus),j+opt.pred_layer_count);
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getReducedModel(optIDX);
                tempMod.init = tempMod.opt;
                xmin = X{1}(1);
                xmax = X{1}(end);
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax);

                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeCenters(opt, X(range), CP(range));
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertPara(tempMod.opt, optIDX);
                SpStatus(j) = 2;
            end
            if obj.options.edgeopt && opt.pred_layer_count < length(Tsidx) % Optimize first and last few layers extra
                for j=2:opt.pred_layer_count
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(1:j)])

                    optIDX = 1:j;

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getReducedModel(optIDX);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeCenters(opt, X(range), CP(range));
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertPara(tempMod.opt, optIDX);
                end
                for j=length(Tsidx)-1:-1:length(Tsidx)-opt.pred_layer_count+1
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(j:length(Tsidx))])

                    optIDX = j:length(Tsidx);

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getReducedModel(optIDX);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeCenters(opt, X(range), CP(range));
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertPara(tempMod.opt, optIDX);
                end
            end
        end
        %% Plotting the current 'CP'-values
        function plotCP(obj)
            T = obj.data.T;
            CP = obj.data.CP;
            figure();
            hold on
            for i=1:length(T)
                plot(CP{i}, T(i)*ones(length(CP{i}),1), 'rx', 'Linewidth', 2, 'Markersize', 16);
            end
            set(gca, 'XDir', 'reverse');
            xlim([obj.data.X{1}(1),obj.data.X{1}(end)]);
            ylim([obj.data.T(1), obj.data.T(end)]);
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variant -1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Make sure that splits are ok and dont stretch or thin the data too much
        function [splits, splitidx, pps] = MakeSplits(obj, t, tempMod, tempind, xmin, xmax)
            numnets = numel(dir(['python/pytorch' '/*.pt'])); % Don't put anything else python-related into this folder!
            ppp = obj.data.CP{t};
            ppp = sort(ppp);
            ppp = ppp(ppp >=xmin & ppp <= xmax); % In the case of var. -2 need to only look in subsystem.
            diffs = diff(ppp);
            splits = zeros(1,length(diffs)+2);
            splits(1) = xmin;
            for i=1:length(diffs) % Split if possible peaks lie > splitppm apart, split in the middle
                if diffs(i) > obj.options.ML.splitppm
                    splits(i+1) = (ppp(i)+ppp(i+1))/2;
                end
            end
            splits(end) = xmax;
            splits = splits(splits ~= 0); % Cuts off unused splits
            splitidx = zeros(size(splits));
            for i=1:length(splits)
                splitidx(i) = find(obj.data.X{1}>=splits(i),1);
            end
            while max(diff(splitidx)) >= 5000 % Then, window is too large and might be compressed too much
                sind = 1;
                while sind < length(splitidx)
                    if splitidx(sind+1)-splitidx(sind) >= 5000
                        ps = [splits(sind), ppp(min(ppp>=splits(sind), ppp<=splits(sind+1)))', splits(sind+1)]; % Current peaks in large window
                        [~, dvind] = max(diff(ps));
                        newsplit = (ps(dvind)+ps(dvind+1))/2;
                        newsplitidx = find(obj.data.X{1}>=newsplit, 1);
                        splits = sort([splits, newsplit]); % Sort new value into old vector
                        splitidx = sort([splitidx, newsplitidx]);
                        sind = sind + 1; % Since a new value was inserted, we advance by two here
                    end
                    sind = sind + 1;
                end
            end
            % What if splits leave too few frequency channels? Then glue together
            pps = zeros(1, length(splits)-1); % Peaks per split
            peakppm = zeros(1, sum(tempMod.init.amt(:,tempind))); % Contains the presumed ppm-values of the multiplet's peaks
            pind = 1;
            for p = 1:size(tempMod.init.amt, 1)
                if tempMod.init.amt(p, tempind) == 1
                    peakppm(pind) = tempMod.init.center(p, tempind);
                    pind = pind + 1;
                else
                    for pp = 1:tempMod.init.amt(p, tempind)
                        peakppm(pind) = tempMod.init.center(p, tempind) + (pp-1)*tempMod.init.dist(p, tempind);
                        pind = pind + 1;
                    end
                end
            end
            for i=1:length(pps)
                pps(i) = numel(find(min(peakppm >= splits(i), peakppm < splits(i+1)))); % Number of peaks in area
            end
            diffs = diff(splitidx);
            for i=length(diffs):-1:2
                if (pps(i-1) + pps(i) <= numnets || pps(i) == 0) && min(diffs(i-1:i))<= 200
                    pps(i-1) = pps(i-1)+pps(i);
                    pps(i) = [];
                    splits(i) = [];
                    splitidx(i) = [];
                    diffs(i-1) = diffs(i-1)+diffs(i);
                    diffs(i) = [];
                end
            end
        end
        %% Big optimization process with ML, variant -1 (only few time layers at a time are optimized)
        function obj = runOptimizationMLPredict(obj, TsStartIdx)
            % TsStartIdx is the Index of Tsidx to start with
            opt = obj.options;
            X = obj.data.X;
            D = obj.data.D;
            ma = obj.data.ma;
            Tsidx = obj.model.Tsidx;

            % SpStatus: 0 ~ unoptimized, 1 ~ upcoming, 2 ~ already optimized
            SpStatus = zeros(size(Tsidx));
            SpStatus(TsStartIdx) = 1;

            % Defining proper pyenv
            pyenv('Version', obj.ppath, 'ExecutionMode','OutOfProcess'); % Correct path must be defined in main.
            % The venv must know torch and numpy.

            % Running initial optimization
            disp(['Run initial: Selected layer # ' num2str(TsStartIdx)])
            optIDX = TsStartIdx;
            tempMod = obj.model.getReducedModel(optIDX);
            tempMod.init = tempMod.opt;
            opttime = tempMod.optimizeModel(opt, X(Tsidx(optIDX)), D(Tsidx(optIDX)), ma(Tsidx(optIDX)), obj.spath, obj.BC); % X(inds) etc. gives cell array, not its contents
            disp(['Computation time: ' num2str(opttime) 's'])
            % Transferring results back to original model
            obj.model.insertPara(tempMod.opt, optIDX);
            SpStatus(optIDX) = 2;

            % Running forward
            for j = TsStartIdx +1 : length(Tsidx)
                disp(['Run forward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = max(TsStartIdx,j-opt.pred_layer_count);
                ide = j; % Indices of layers to be accounted for while optimizing in this window
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % First, linearly predict parameters as before (to deal with crossings)
                tempMod = obj.model.getReducedModel(optIDX);
                tempMod.init = tempMod.opt;
                xmin = X{1}(1);
                xmax = X{1}(end);
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax);
                tempind = find(SpStatus(optIDX)==1);
                % Predicting parameters with Machine Learning
                % Split the given spectrum into windows with acceptable sizes and peak numbers
                [splits, splitidx, pps] = obj.MakeSplits(Tsidx(j), tempMod, tempind, xmin, xmax);
                numnets = numel(dir(['python/pytorch' '/*.pt'])); % Counts number of nets, resp. state dicts in folder! Careful not to put unused nets in folder!
                nextspec = D{Tsidx(j)};
                paracell = cell(1,length(splits)-1);
                longspec = zeros(1, (length(splits)-1)*1000);
                heightfacs = ones(1, length(splits)-1);
                for i=1:length(splits)-1
                    specwindow = nextspec(splitidx(i):splitidx(i+1));
                    Xwindow = X{1}(splitidx(i):splitidx(i+1));
                    xnew = linspace(Xwindow(1), Xwindow(end), 1000);
                    Dnew = pchip(Xwindow, specwindow, xnew);
                    if max(Dnew) < 0.75 && max(Dnew) > opt.L1.epsi % Norms the splits, if very small peaks are present only
                        heightfacs(i) = 1/max(Dnew)*0.95; % Nets prefer peaks < 1
                        Dnew = Dnew*heightfacs(i);
                    end
                    longspec(1000*(i-1)+1:1000*i) = Dnew; % Longspec might not produce a sensible spectrum, it is used to store spectra only
                end
                pydict = pyrunfile("python/predict_window3.py", "out_d", longwin=longspec); % Only runs once, more efficient.
                parastruct = struct(pydict); % Params are behind 's<i>p<k>, where 'i' = no. split, 'k' = no. peak predicted
                for i=1:length(splits)-1
                    errs = zeros(1, numnets+1);
                    params = cell(1, numnets);
                    xerr = linspace(0,10,1000);
                    Dnew = longspec(1000*(i-1)+1:1000*i); % Might be scaled additionally
                    for k = 1:numnets
                        name = 's' + string(i) + 'p' + string(k);
                        pas = parastruct.(name); % Dynamic field name, no eval() needed 
                        params{k} = double(pas); % Centers, heights, hws, lambdas
                        Derr = zeros(size(xerr));
                        for np= 1:k
                            Derr = Derr + gl(xerr, 1/2*params{k}(2*k+np), 10*params{k}(np), params{k}(k+np), 1/5*params{k}(3*k+np));
                        end
                        errs(k) = norm(Dnew-Derr);
                    end
                    % Find out how the extrapolated prediction compares
                    pind = find(tempMod.init.center(:,tempind) <= splits(i+1) & tempMod.init.center(:, tempind) >= splits(i)); % Finds indices of peaks in current range
                    Dex = zeros(size(Derr));
                    for np = 1:length(pind)
                        camt = tempMod.init.amt(pind(np),tempind);
                        ccent = 10*(tempMod.init.center(pind(np), tempind)-splits(i))/(splits(i+1)-splits(i));
                        cdist = max(-1, 10*tempMod.init.dist(pind(np), tempind)/(splits(i+1)-splits(i)));
                        chw  = 10*tempMod.init.hw(pind(np), tempind)/(splits(i+1)-splits(i));
                        cheight = tempMod.init.height(pind(np), tempind)*heightfacs(i);
                        clamb = tempMod.init.lambda(pind(np), tempind);
                        Dex = Dex + glmult(xerr, camt, ccent, cdist, chw, cheight, clamb, obj.BC);
                    end
                    exprederr = norm(Dnew-Dex);
                    errs(end) = min(norm(Dnew), exprederr);
                    if pps(i) == 0 || pps(i) > numnets
                        errs(end) = errs(end)*opt.ML.peakadv;
                    else
                        errs(pps(i)) = errs(pps(i))*opt.ML.peakadv; % Else let modif. best fit decide
                    end
                    [~, numpeak] = min(errs);
                    if numpeak ~= length(errs) % Then peaks are present, else there aren't or every prediction is worse than not predicting
                        left = splits(i);
                        right = splits(i+1);
                        sz = right-left;
                        temppara = params{numpeak};
                        temppara(1:numpeak) = temppara(1:numpeak)*sz + left; % Scale centers
                        temppara(numpeak+1:2*numpeak) = 1/heightfacs(i)*temppara(numpeak+1:2*numpeak); % Scale heights back
                        temppara(2*numpeak+1:3*numpeak) = sz/20*temppara(2*numpeak+1:3*numpeak); % Scale hws (1/2*1/10*sz)
                        temppara(3*numpeak+1:end) = 1/5*temppara(3*numpeak+1:end); % Scale lambdas
                        paracell{i} = temppara;
                    end
                end

                % After paracell has been obtained sort obtained parameters and peaks into next layer of tempmodel
                peakppm = zeros(1, sum(tempMod.init.amt(:,tempind))); % Assumed ppm-values of multiplet's subpeaks
                multidx = peakppm; % Memory which peak belonged to which multiplet
                mnum = peakppm; % Memory where in multiplet current peak stands
                pind = 1;
                for p = 1:size(tempMod.init.amt, 1)
                    if tempMod.init.amt(p, tempind) == 1
                        peakppm(pind) = tempMod.init.center(p, tempind);
                        multidx(pind) = p;
                        mnum(pind) = 1;
                        pind = pind + 1;
                    else
                        for pp = 1:tempMod.init.amt(p, tempind)
                            peakppm(pind) = tempMod.init.center(p, tempind) + (pp-1)*tempMod.init.dist(p, tempind);
                            multidx(pind) = p;
                            mnum(pind) = pp;
                            pind = pind + 1;
                        end
                    end
                end
                tempstruct = tempMod.init; % Allowed to overwrite things in tempstruct
                fns = fieldnames(tempstruct);
                for n = [2,4,5,6] % Only need to overwrite certain fields in the struct
                    tempstruct.(fns{n}) = zeros(size(tempstruct.(fns{n})));
                end
                tempstruct.prednum = zeros(size(tempstruct.center,1),1);
                for ind = 1:length(paracell)
                    if ~isempty(paracell{ind})
                        logarr = min(peakppm>=splits(ind),peakppm<=splits(ind+1)); % Minimum represents logical AND
                        splitc = peakppm(logarr); 
                        splitinds = multidx(logarr);
                        smnum = mnum(logarr);
                        paramat = reshape(paracell{ind}, [], 4)'; % [centers; heights; hws; lambdas]
                        [~, order] = sort(paramat(1,:));
                        paramat(:, :) = paramat(:, order); % Sorts matrix by first row (should be sorted, but with ML you never know)
                        im = length(splitc);
                        jm = size(paramat,2);
                        distmat = zeros(im, jm); % We calculate the distances of all predicted peaks to the extrapolated peaks and sort the nearest first
                        for iim = 1:im
                            for jjm = 1:jm
                                distmat(iim, jjm) = abs(splitc(iim)-paramat(1,jjm));
                            end
                        end
                        while im > 0 && jm > 0
                            [Ai, Aj] = find(distmat == min(distmat(:)),1);
                            distmat(Ai, :) = ones(size(distmat(Ai, :)))*inf;
                            distmat(:, Aj) = ones(size(distmat(:, Aj)))*inf;
                            im = im-1;
                            jm = jm-1;
                            nump = smnum(Ai);
                            Ai = splitinds(Ai);
                            % We will take the average over all predictions of the corresponding multiplet parameters
                            if nump > 1
                                tempstruct.center(Ai, tempind) = tempstruct.center(Ai, tempind) + paramat(1, Aj) - tempstruct.dist(Ai, tempind)*(nump-1); % Adjust centers to the fist peak of the multiplet
                                tempstruct.height(Ai, tempind) = tempstruct.height(Ai, tempind) + paramat(2, Aj)/obj.BC(tempstruct.amt(Ai, tempind), nump); % Adjust heights
                            else
                                tempstruct.center(Ai, tempind) = tempstruct.center(Ai, tempind) + paramat(1, Aj);
                                tempstruct.height(Ai, tempind) = tempstruct.height(Ai, tempind) + paramat(2, Aj);
                            end
                            tempstruct.hw(Ai, tempind) = tempstruct.hw(Ai, tempind) + paramat(3, Aj);
                            tempstruct.lambda(Ai, tempind) = tempstruct.lambda(Ai, tempind) + paramat(4, Aj);
                            tempstruct.prednum(Ai) = tempstruct.prednum(Ai)+1;
                        end
                    end
                end
                % Prevent division by 0
                lvec = tempstruct.prednum > 0; % No changes in lvec, no need to change the init
                for n = [2 4 5 6] % Averaging and initializing tempMod
                    tempMod.init.(fns{n})(lvec, tempind) = tempstruct.(fns{n})(lvec, tempind)./tempstruct.prednum(lvec);
                end
                % More plots
                if opt.info == 2.2
                    disp(['Plot initial prediction of layer # ' num2str(Tsidx(j))])
                    d = D{Tsidx(j)};
                    x = X{Tsidx(j)};
                    MLmodel = zeros(size(x));
                    for p = 1:size(tempMod.init.center, 1)
                        MLmodel = MLmodel + glmult(x, tempMod.init.amt(p, tempind), tempMod.init.center(p, tempind), tempMod.init.dist(p, tempind), tempMod.init.hw(p, tempind), tempMod.init.height(p, tempind), tempMod.init.lambda(p, tempind), obj.BC);
                    end
                    figure(7);
                    set(gcf, 'Position', get(0, 'Screensize'));
                    if opt.epsinits
                        clf
                        hold on
                        plot(x,d,'k', 'LineWidth', 2)
                        plot(x,MLmodel, 'b', 'LineWidth', 2)
                    else
                        plot(x, d, 'k', x, MLmodel, 'b');
                    end
                    set(gca, 'XDir', 'reverse')
                    set(gcf, 'renderer', 'painters');
                    title('t1')
                    xlabel('x1')
                    ylabel('y1')
                    set(gca, 'Fontsize', 24)
                    if opt.epsinits
                        ymax = max(max(d),max(MLmodel));
                        ylim([0,ymax])
                        set(gca, 'Fontsize', 48)
                        path = ['Bilder/Ints/abc125', num2str(Tsidx(j))];
                        pause(2)
                        print('-depsc', '-loose', path);
                        hold off
                    end
                    wm = wasserstein(x, MLmodel, d)*100/(x(end)-x(1)); % Normalized Wasserstein metric for one spectrum
                    frobm = norm(MLmodel - d, 'fro')*100/length(x); % MSE
                    disp(['Wasserstein-Metric of: ' num2str(wm) ' %'])
                    disp(['MSE of:' num2str(frobm) ' %'])
                end
                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertPara(tempMod.opt, optIDX);
                SpStatus(j) = 2;
            end
            % Running backwards
            for j = TsStartIdx -1:-1:1
                disp(['Run backward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = j;
                ide = min(length(SpStatus),j+opt.pred_layer_count);
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % First, linearly predict parameters as before (to deal with crossings)
                tempMod = obj.model.getReducedModel(optIDX);
                tempMod.init = tempMod.opt;
                xmin = X{1}(1);
                xmax = X{1}(end);
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax);
                tempind = find(SpStatus(optIDX)==1);
                % Predicting parameters with Machine Learning
                % Split the given spectrum into windows with acceptable sizes and peak numbers
                [splits, splitidx, pps] = obj.MakeSplits(Tsidx(j), tempMod, tempind, xmin, xmax);
                numnets = numel(dir(['python/pytorch' '/*.pt'])); % Counts number of nets, resp. state dicts in folder! Careful not to put unused nets in folder!
                nextspec = D{Tsidx(j)};
                paracell = cell(1,length(splits)-1);
                longspec = zeros(1, (length(splits)-1)*1000);
                heightfacs = ones(1, length(splits)-1);
                for i=1:length(splits)-1
                    specwindow = nextspec(splitidx(i):splitidx(i+1));
                    Xwindow = X{1}(splitidx(i):splitidx(i+1));
                    xnew = linspace(Xwindow(1), Xwindow(end), 1000);
                    Dnew = pchip(Xwindow, specwindow, xnew);
                    if max(Dnew) < 0.75 && max(Dnew) > opt.L1.epsi % Norms the splits, if very small peaks are present only
                        heightfacs(i) = 1/max(Dnew)*0.95; % Nets prefer peaks < 1
                        Dnew = Dnew*heightfacs(i);
                    end
                    longspec(1000*(i-1)+1:1000*i) = Dnew; % Longspec might not produce a sensible spectrum, it is used to store spectra only
                end
                pydict = pyrunfile("python/predict_window3.py", "out_d", longwin=longspec); % Only runs once
                parastruct = struct(pydict); % Params are behind "s<i>p<k>", 'i' = no. split, 'k' = no. peak predicted
                for i=1:length(splits)-1
                    errs = zeros(1, numnets+1);
                    params = cell(1, numnets);
                    xerr = linspace(0,10,1000);
                    Dnew = longspec(1000*(i-1)+1:1000*i); % Might be scaled additionally
                    for k = 1:numnets
                        name = 's' + string(i) + 'p' + string(k);
                        pas = parastruct.(name); % Dynamic field name, no eval() needed
                        params{k} = double(pas); % Centers, heights, hws, lambdas
                        Derr = zeros(size(xerr));
                        for np= 1:k
                            Derr = Derr + gl(xerr, 1/2*params{k}(2*k+np), 10*params{k}(np), params{k}(k+np), 1/5*params{k}(3*k+np));
                        end
                        errs(k) = norm(Dnew-Derr);
                    end
                    % Find out how the extrapolated prediction compares
                    pind = find(tempMod.init.center(:,tempind) <= splits(i+1) & tempMod.init.center(:, tempind) >= splits(i)); % Finds indices of peaks in current range
                    Dex = zeros(size(Derr));
                    for np = 1:length(pind)
                        camt = tempMod.init.amt(pind(np),tempind);
                        ccent = 10*(tempMod.init.center(pind(np), tempind)-splits(i))/(splits(i+1)-splits(i));
                        cdist = max(-1, 10*tempMod.init.dist(pind(np), tempind)/(splits(i+1)-splits(i)));
                        chw  = 10*tempMod.init.hw(pind(np), tempind)/(splits(i+1)-splits(i));
                        cheight = tempMod.init.height(pind(np), tempind)*heightfacs(i);
                        clamb = tempMod.init.lambda(pind(np), tempind);
                        Dex = Dex + glmult(xerr, camt, ccent, cdist, chw, cheight, clamb, obj.BC);
                    end
                    exprederr = norm(Dnew-Dex);
                    errs(end) = min(norm(Dnew), exprederr);
                    if pps(i) == 0 || pps(i) > numnets
                        errs(end) = errs(end)*opt.ML.peakadv;
                    else
                        errs(pps(i)) = errs(pps(i))*opt.ML.peakadv; % Else let modif. best fit decide
                    end
                    [~, numpeak] = min(errs);
                    if numpeak ~= length(errs) % Then peaks are present, else there aren't or every prediction is worse than not predicting.
                        left = splits(i);
                        right = splits(i+1);
                        sz = right-left;
                        temppara = params{numpeak};
                        temppara(1:numpeak) = temppara(1:numpeak)*sz + left; % Scale centers
                        temppara(numpeak+1:2*numpeak) = 1/heightfacs(i)*temppara(numpeak+1:2*numpeak); % Scale heights back
                        temppara(2*numpeak+1:3*numpeak) = sz/20*temppara(2*numpeak+1:3*numpeak); % Scale hws (1/2*1/10*sz)
                        temppara(3*numpeak+1:end) = 1/5*temppara(3*numpeak+1:end); % Scale lambdas
                        paracell{i} = temppara;
                    end
                end

                % After paracell has been obtained need to sort obtained parameters and peaks into next layer of tempmodel
                peakppm = zeros(1, sum(tempMod.init.amt(:,tempind))); % Assumed ppm-values of multiplet's subpeaks
                multidx = peakppm; % Memory which peak belonged to which multiplet
                mnum = peakppm; % Memory where in multiplet current peak stands
                pind = 1;
                for p = 1:size(tempMod.init.amt, 1)
                    if tempMod.init.amt(p, tempind) == 1
                        peakppm(pind) = tempMod.init.center(p, tempind);
                        multidx(pind) = p;
                        mnum(pind) = 1;
                        pind = pind + 1;
                    else
                        for pp = 1:tempMod.init.amt(p, tempind)
                            peakppm(pind) = tempMod.init.center(p, tempind) + (pp-1)*tempMod.init.dist(p, tempind);
                            multidx(pind) = p;
                            mnum(pind) = pp;
                            pind = pind + 1;
                        end
                    end
                end
                tempstruct = tempMod.init; % Allowed to overwrite things in tempstruct
                fns = fieldnames(tempstruct);
                for n = [2,4,5,6] % Only need to overwrite certain fields in the struct
                    tempstruct.(fns{n}) = zeros(size(tempstruct.(fns{n})));
                end
                tempstruct.prednum = zeros(size(tempstruct.center,1),1);
                for ind = 1:length(paracell)
                    if ~isempty(paracell{ind})
                        logarr = min(peakppm>=splits(ind),peakppm<=splits(ind+1)); % Minimum represents logical AND
                        splitc = peakppm(logarr); 
                        splitinds = multidx(logarr);
                        smnum = mnum(logarr);
                        paramat = reshape(paracell{ind}, [], 4)'; % [centers; heights; hws; lambdas]
                        [~, order] = sort(paramat(1,:));
                        paramat(:, :) = paramat(:, order); % Sorts matrix by first row (should be sorted, but with ML you never know)
                        im = length(splitc);
                        jm = size(paramat,2);
                        distmat = zeros(im, jm); % We calculate the distances of all predicted peaks to the extrapolated peaks and sort the nearest first
                        for iim = 1:im
                            for jjm = 1:jm
                                distmat(iim, jjm) = abs(splitc(iim)-paramat(1,jjm));
                            end
                        end
                        while im > 0 && jm > 0
                            [Ai, Aj] = find(distmat == min(distmat(:)),1);
                            distmat(Ai, :) = ones(size(distmat(Ai, :)))*inf;
                            distmat(:, Aj) = ones(size(distmat(:, Aj)))*inf;
                            im = im-1;
                            jm = jm-1;
                            nump = smnum(Ai);
                            Ai = splitinds(Ai);
                            % We take the average over all predictions of the corresponding multiplet parameters
                            if nump > 1
                                tempstruct.center(Ai, tempind) = tempstruct.center(Ai, tempind) + paramat(1, Aj) - tempstruct.dist(Ai, tempind)*(nump-1); % Adjust centers to the fist peak of the multiplet
                                tempstruct.height(Ai, tempind) = tempstruct.height(Ai, tempind) + paramat(2, Aj)/obj.BC(tempstruct.amt(Ai, tempind), nump); % Adjust heights
                            else
                                tempstruct.center(Ai, tempind) = tempstruct.center(Ai, tempind) + paramat(1, Aj);
                                tempstruct.height(Ai, tempind) = tempstruct.height(Ai, tempind) + paramat(2, Aj);
                            end
                            tempstruct.hw(Ai, tempind) = tempstruct.hw(Ai, tempind) + paramat(3, Aj);
                            tempstruct.lambda(Ai, tempind) = tempstruct.lambda(Ai, tempind) + paramat(4, Aj);
                            tempstruct.prednum(Ai) = tempstruct.prednum(Ai)+1;
                        end
                    end
                end
                % Prevent division by 0.
                lvec = tempstruct.prednum > 0; % No changes in lvec, no need to change the init
                for n = [2 4 5 6] % Averaging and initializing tempMod
                    tempMod.init.(fns{n})(lvec, tempind) = tempstruct.(fns{n})(lvec, tempind)./tempstruct.prednum(lvec);
                end
                % Even more plots
                if opt.info == 2.2
                    disp(['Plot initial prediction of layer # ' num2str(Tsidx(j))])
                    d = D{Tsidx(j)};
                    x = X{Tsidx(j)};
                    MLmodel = zeros(size(x));
                    for p = 1:size(tempMod.init.center, 1)
                        MLmodel = MLmodel + glmult(x, tempMod.init.amt(p, tempind), tempMod.init.center(p, tempind), tempMod.init.dist(p, tempind), tempMod.init.hw(p, tempind), tempMod.init.height(p, tempind), tempMod.init.lambda(p, tempind), obj.BC);
                    end
                    figure(7);
                    set(gcf, 'Position', get(0, 'Screensize'));
                    if opt.epsinits
                        clf
                        hold on
                        plot(x,d,'k', 'LineWidth', 2)
                        plot(x,MLmodel, 'b', 'LineWidth', 2)
                    else
                        plot(x, d, 'k', x, MLmodel, 'b');
                    end
                    set(gca, 'XDir', 'reverse')
                    set(gcf, 'renderer', 'painters');
                    title('t1')
                    xlabel('x1')
                    ylabel('y1')
                    set(gca, 'Fontsize', 24)
                    if opt.epsinits
                        ymax = max(max(d),max(MLmodel));
                        ylim([0,ymax])
                        set(gca, 'Fontsize', 48)
                        path = ['Bilder/Ints/abc126', num2str(Tsidx(j))];
                        pause(2)
                        print('-depsc', '-loose', path);
                        hold off
                    end
                    wm = wasserstein(x, MLmodel, d)*100/(x(end)-x(1)); % Normalized Wasserstein metric for one spectrum
                    frobm = norm(MLmodel - d, 'fro')*100/length(x); % MSE
                    disp(['Wasserstein-Metric of: ' num2str(wm) ' %'])
                    disp(['MSE of:' num2str(frobm) ' %'])
                end
                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertPara(tempMod.opt, optIDX);
                SpStatus(j) = 2;
            end
            if obj.options.edgeopt && opt.pred_layer_count < length(Tsidx) % Optimize first and last few layers extra
                for j=2:opt.pred_layer_count
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(1:j)])

                    optIDX = 1:j;

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getReducedModel(optIDX);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertPara(tempMod.opt, optIDX);
                end
                for j=length(Tsidx)-1:-1:length(Tsidx)-opt.pred_layer_count+1
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(j:length(Tsidx))])

                    optIDX = j:length(Tsidx);

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getReducedModel(optIDX);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertPara(tempMod.opt, optIDX);
                end
            end
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variant -2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Finding boudaries to split the data set using heuristic Peak detection only!
        function [splits, multidx] = findOptboundsML(obj)
            X = obj.data.X{1};
            splits = X(1);
            CP = obj.data.CP; % Copy, for rearranging later
            s = 0;
            for t=1:length(CP)
                CP{t} = sort(CP{t});
                s = s + length(CP{t});
            end
            CPvec = zeros(1, s);
            ind = 1;
            for t = 1:length(CP)
                CPvec(ind:ind+length(CP{t}) - 1) = CP{t};
                ind = ind + length(CP{t});
            end
            CPvec = sort(CPvec);
            CPdiff = diff(CPvec);
            for d = 1:length(CPdiff)
                if CPdiff(d) > obj.options.L1.splitdist
                    splitpoint = (CPvec(d+1)+CPvec(d))/2; % Choose midpoint to split
                    splitidx = find(X >= splitpoint, 1);
                    splits = [splits, X(splitidx)]; %#ok<AGROW>
                end
            end
            splits = [splits, X(end)];
            % Sort multiplets into splits by using initial centers
            cinit = obj.model.init.center(:,1);
            multidx = cell(1, length(splits)-1);
            for s = 1:length(multidx)
                multidx{s} = find(cinit >= splits(s) & cinit <= splits(s+1));
            end
        end
        %% Optimizing all variables occurring in certain bounds
        function obj = runBoundedOptML(obj, TsStartIdx, lb, ub, multidx)
            % TsStartIdx is the Index of Tsidx to start with
            opt = obj.options;
            X = obj.data.X;
            D = obj.data.D;
            xmin = lb;
            xmax = ub;
            Tsidx = obj.model.Tsidx;
            ma = obj.data.ma;

            % Limiting X and D
            lbind = find(X{1}>= lb, 1, 'first');
            ubind = find(X{1}<= ub, 1, 'last');

            for t = 1:length(X)
                X{t} = X{t}(lbind:ubind); % Not permanent on data.X :)
                D{t} = D{t}(lbind:ubind);
                ma(t) = max(abs(D{t}));
            end

            % SpStatus: 0 ~ unoptimized, 1 ~ upcoming, 2 ~ already optimized
            SpStatus = zeros(size(Tsidx));
            SpStatus(TsStartIdx) = 1;

            % Defining proper pyenv
            pyenv('Version',obj.ppath, 'ExecutionMode','OutOfProcess'); % Include proper path defined in main.
            % Venv must know about torch and numpy

            % Running initial optimization
            disp(['Run initial: Selected layer # ' num2str(TsStartIdx)])
            optIDX = TsStartIdx;
            tempMod = obj.model.getBoundedModel(optIDX, multidx); % Needs creating
            tempMod.init = tempMod.opt; % Needs corrections
            tempMod.lower.center = max(tempMod.lower.center, lb*ones(size(tempMod.lower.center))); % If somehow stricter limits exist, keep those
            tempMod.upper.center = min(tempMod.upper.center, ub*ones(size(tempMod.upper.center)));
            opttime = tempMod.optimizeModel(opt, X(Tsidx(optIDX)), D(Tsidx(optIDX)), ma(Tsidx(optIDX)), obj.spath, obj.BC); % 'X(inds)' etc. gives cell array, not its contents
            disp(['Computation time: ' num2str(opttime) 's'])
            % Transferring results back to original model
            obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx); % Needs adapting
            SpStatus(optIDX) = 2;

            % Running forward
            for j = TsStartIdx +1 : length(Tsidx)
                disp(['Run forward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = max(TsStartIdx,j-opt.pred_layer_count);
                ide = j;
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getBoundedModel(optIDX, multidx);
                tempMod.init = tempMod.opt;
                tempMod.lower.center = lb*ones(size(tempMod.lower.center));
                tempMod.upper.center = ub*ones(size(tempMod.upper.center));

                tempMod.interpPara(SpStatus(optIDX), xmin, xmax)
                tempind = find(SpStatus(optIDX)==1);
                % Predicting parameters with Machine Learning
                % Split the given spectrum into windows with acceptable sizes and peak numbers
                [splits, splitidx, pps] = obj.MakeSplits(Tsidx(j), tempMod, tempind, xmin, xmax);
                splitidx = splitidx - splitidx(1) + 1; % Calculate indices of current X, not global X.
                numnets = numel(dir(['python/pytorch' '/*.pt'])); % Counts number of nets, resp. state dicts in folder! Careful not to put unused nets in folder!
                nextspec = D{Tsidx(j)};
                if length(nextspec) + 1 == splitidx(end)
                    splitidx(end) = length(nextspec);
                end
                paracell = cell(1,length(splits)-1);
                longspec = zeros(1, (length(splits)-1)*1000);
                heightfacs = ones(1, length(splits)-1);
                for i=1:length(splits)-1
                    specwindow = nextspec(splitidx(i):splitidx(i+1));
                    Xwindow = X{1}(splitidx(i):splitidx(i+1));
                    xnew = linspace(Xwindow(1), Xwindow(end), 1000);
                    Dnew = pchip(Xwindow, specwindow, xnew);
                    if max(Dnew) < 0.75 && max(Dnew) > opt.L1.epsi % Normalizes the splits, if small peaks are present only
                        heightfacs(i) = 1/max(Dnew)*0.95; % Nets prefer peaks < 1
                        Dnew = Dnew*heightfacs(i);
                    end
                    longspec(1000*(i-1)+1:1000*i) = Dnew; % Longspec might not produce a sensible spectrum, it is used to store spectra only
                end
                pydict = pyrunfile("python/predict_window3.py", "out_d", longwin=longspec); % Only runs once
                parastruct = struct(pydict); % Params are behind "s<i>p<k>", 'i' = no. split, 'k' = no. peaks predicted
                for i=1:length(splits)-1
                    errs = zeros(1, numnets+1);
                    params = cell(1, numnets);
                    xerr = linspace(0,10,1000);
                    Dnew = longspec(1000*(i-1)+1:1000*i); % Might be scaled additionally
                    for k = 1:numnets
                        name = 's' + string(i) + 'p' + string(k);
                        pas = parastruct.(name); % Dynamic field name, no eval() needed
                        params{k} = double(pas); % Centers, heights, hws, lambdas
                        Derr = zeros(size(xerr));
                        for np= 1:k
                            Derr = Derr + gl(xerr, 1/2*params{k}(2*k+np), 10*params{k}(np), params{k}(k+np), 1/5*params{k}(3*k+np));
                        end
                        errs(k) = norm(Dnew-Derr);
                    end
                    % Find out how the extrapolated prediction compares
                    pind = find(tempMod.init.center(:,tempind) <= splits(i+1) & tempMod.init.center(:, tempind) >= splits(i)); % finds indices of peaks in current range
                    Dex = zeros(size(Derr));
                    for np = 1:length(pind)
                        camt = tempMod.init.amt(pind(np),tempind);
                        ccent = 10*(tempMod.init.center(pind(np), tempind)-splits(i))/(splits(i+1)-splits(i));
                        cdist = max(-1, 10*tempMod.init.dist(pind(np), tempind)/(splits(i+1)-splits(i)));
                        chw  = 10*tempMod.init.hw(pind(np), tempind)/(splits(i+1)-splits(i));
                        cheight = tempMod.init.height(pind(np), tempind)*heightfacs(i);
                        clamb = tempMod.init.lambda(pind(np), tempind);
                        Dex = Dex + glmult(xerr, camt, ccent, cdist, chw, cheight, clamb, obj.BC);
                    end
                    exprederr = norm(Dnew-Dex);
                    errs(end) = min(norm(Dnew), exprederr);
                    if pps(i) == 0 || pps(i) > numnets
                        errs(end) = errs(end)*opt.ML.peakadv;
                    else
                        errs(pps(i)) = errs(pps(i))*opt.ML.peakadv; % Else let modif. best fit decide
                    end
                    [~, numpeak] = min(errs);
                    if numpeak ~= length(errs) % Then peaks are present, else they aren't, or every prediction is worse than not predicting.
                        left = splits(i);
                        right = splits(i+1);
                        sz = right-left;
                        temppara = params{numpeak};
                        temppara(1:numpeak) = temppara(1:numpeak)*sz + left; % Scale centers
                        temppara(numpeak+1:2*numpeak) = 1/heightfacs(i)*temppara(numpeak+1:2*numpeak); % Scale heights back
                        temppara(2*numpeak+1:3*numpeak) = sz/20*temppara(2*numpeak+1:3*numpeak); % Scale hws (1/2*1/10*sz)
                        temppara(3*numpeak+1:end) = 1/5*temppara(3*numpeak+1:end); % Scale lambdas
                        paracell{i} = temppara;
                    end
                end

                % After paracell has been obtained need to sort obtained parameters and peaks into next layer of tempmodel
                peakppm = zeros(1, sum(tempMod.init.amt(:,tempind))); % Assumed ppm-values of multiplet's subpeaks
                multipidx = peakppm; % Memory which peak belonged to which multiplet
                mnum = peakppm; % Memory where in multiplet current peak stands
                pind = 1;
                for p = 1:size(tempMod.init.amt, 1)
                    if tempMod.init.amt(p, tempind) == 1
                        peakppm(pind) = tempMod.init.center(p, tempind);
                        multipidx(pind) = p;
                        mnum(pind) = 1;
                        pind = pind + 1;
                    else
                        for pp = 1:tempMod.init.amt(p, tempind)
                            peakppm(pind) = tempMod.init.center(p, tempind) + (pp-1)*tempMod.init.dist(p, tempind);
                            multipidx(pind) = p;
                            mnum(pind) = pp;
                            pind = pind + 1;
                        end
                    end
                end
                tempstruct = tempMod.init; % Allowed to overwrite things in tempstruct
                fns = fieldnames(tempstruct);
                for n = [2,4,5,6] % Only need to overwrite certain fields in the struct
                    tempstruct.(fns{n}) = zeros(size(tempstruct.(fns{n})));
                end
                tempstruct.prednum = zeros(size(tempstruct.center,1),1);
                for ind = 1:length(paracell)
                    if ~isempty(paracell{ind})
                        logarr = min(peakppm>=splits(ind),peakppm<=splits(ind+1)); % Minimum represents logical AND
                        splitc = peakppm(logarr); 
                        splitinds = multipidx(logarr);
                        smnum = mnum(logarr);
                        paramat = reshape(paracell{ind}, [], 4)'; % [centers; heights; hws; lambdas]
                        [~, order] = sort(paramat(1,:));
                        paramat(:, :) = paramat(:, order); % Sorts matrix by first row (should be sorted, but with ML you never know)
                        im = length(splitc);
                        jm = size(paramat,2);
                        distmat = zeros(im, jm); % We calculate the distances of all predicted peaks to the extrapolated peaks and sort the nearest first
                        for iim = 1:im
                            for jjm = 1:jm
                                distmat(iim, jjm) = abs(splitc(iim)-paramat(1,jjm));
                            end
                        end
                        while im > 0 && jm > 0
                            [Ai, Aj] = find(distmat == min(distmat(:)),1);
                            distmat(Ai, :) = ones(size(distmat(Ai, :)))*inf;
                            distmat(:, Aj) = ones(size(distmat(:, Aj)))*inf;
                            im = im-1;
                            jm = jm-1;
                            nump = smnum(Ai);
                            Ai = splitinds(Ai);
                            % We take the average over all predictions of the corresponding multiplet parameters
                            if nump > 1
                                tempstruct.center(Ai, tempind) = tempstruct.center(Ai, tempind) + paramat(1, Aj) - tempstruct.dist(Ai, tempind)*(nump-1); % adjust centers to the fist peak of the multiplet
                                tempstruct.height(Ai, tempind) = tempstruct.height(Ai, tempind) + paramat(2, Aj)/obj.BC(tempstruct.amt(Ai, tempind), nump); % adjust heights
                            else
                                tempstruct.center(Ai, tempind) = tempstruct.center(Ai, tempind) + paramat(1, Aj);
                                tempstruct.height(Ai, tempind) = tempstruct.height(Ai, tempind) + paramat(2, Aj);
                            end
                            tempstruct.hw(Ai, tempind) = tempstruct.hw(Ai, tempind) + paramat(3, Aj);
                            tempstruct.lambda(Ai, tempind) = tempstruct.lambda(Ai, tempind) + paramat(4, Aj);
                            tempstruct.prednum(Ai) = tempstruct.prednum(Ai)+1;
                        end
                    end
                end
                % Prevents division by 0.
                lvec = tempstruct.prednum > 0; % No changes in lvec, no need to change the init
                for n = [2 4 5 6] % Averaging and initializing tempMod
                    tempMod.init.(fns{n})(lvec, tempind) = tempstruct.(fns{n})(lvec, tempind)./tempstruct.prednum(lvec);
                end
                % Further plots
                if opt.info == 2.2
                    disp(['Plot initial prediction of layer # ' num2str(Tsidx(j))])
                    d = D{Tsidx(j)};
                    x = X{Tsidx(j)};
                    MLmodel = zeros(size(x));
                    for p = 1:size(tempMod.init.center, 1)
                        MLmodel = MLmodel + glmult(x, tempMod.init.amt(p, tempind), tempMod.init.center(p, tempind), tempMod.init.dist(p, tempind), tempMod.init.hw(p, tempind), tempMod.init.height(p, tempind), tempMod.init.lambda(p, tempind), obj.BC);
                    end
                    figure(7);
                    set(gcf, 'Position', get(0, 'Screensize'));
                    if opt.epsinits
                        clf
                        hold on
                        plot(x,d,'k', 'LineWidth', 2)
                        plot(x,MLmodel, 'b', 'LineWidth', 2)
                    else
                        plot(x, d, 'k', x, MLmodel, 'b');
                    end
                    set(gca, 'XDir', 'reverse')
                    set(gcf, 'renderer', 'painters');
                    title('t1')
                    xlabel('x1')
                    xlim([xmin, xmax])
                    ylabel('y1')
                    set(gca, 'Fontsize', 24)
                    if opt.epsinits
                        ymax = max(max(d),max(MLmodel));
                        ylim([0,ymax])
                        set(gca, 'Fontsize', 48)
                        path = ['Bilder/Inits/abc127', num2str(Tsidx(j))];
                        pause(2)
                        print('-depsc', '-loose', path);
                        hold off
                    end
                    wm = wasserstein(x, MLmodel, d)*100/(x(end)-x(1)); % Normalized Wasserstein metric for one spectrum
                    frobm = norm(MLmodel - d, 'fro')*100/length(x); % MSE
                    disp(['Wasserstein-Metric of: ' num2str(wm) ' %'])
                    disp(['MSE of:' num2str(frobm) ' %'])
                end
                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                SpStatus(j) = 2;
            end

            % Running backwards
            for j = TsStartIdx-1:-1:1
                disp(['Run backward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = j;
                ide = min(length(SpStatus),j+opt.pred_layer_count);
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getBoundedModel(optIDX, multidx);
                tempMod.init = tempMod.opt;
                tempMod.lower.center = lb*ones(size(tempMod.lower.center));
                tempMod.upper.center = ub*ones(size(tempMod.upper.center));
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax);
                tempind = find(SpStatus(optIDX)==1);

                % Predicting parameters with Machine Learning
                % Split the given spectrum into windows with acceptable sizes and peak numbers
                [splits, splitidx, pps] = obj.MakeSplits(Tsidx(j), tempMod, tempind, xmin, xmax);
                splitidx = splitidx - splitidx(1) + 1; % Indices of current X, not global X!
                numnets = numel(dir(['python/pytorch' '/*.pt'])); % Counts number of nets, resp. state dicts in folder! Careful not to put unused nets in folder!
                nextspec = D{Tsidx(j)};
                if length(nextspec) + 1 == splitidx(end)
                    splitidx(end) = length(nextspec);
                end
                paracell = cell(1,length(splits)-1);
                longspec = zeros(1, (length(splits)-1)*1000);
                heightfacs = ones(1, length(splits)-1);
                for i=1:length(splits)-1
                    specwindow = nextspec(splitidx(i):splitidx(i+1));
                    Xwindow = X{1}(splitidx(i):splitidx(i+1));
                    xnew = linspace(Xwindow(1), Xwindow(end), 1000);
                    Dnew = pchip(Xwindow, specwindow, xnew);
                    if max(Dnew) < 0.75 && max(Dnew) > opt.L1.epsi % Norms the splits, if very small peaks are present only
                        heightfacs(i) = 1/max(Dnew)*0.95; % Nets prefer peaks < 1
                        Dnew = Dnew*heightfacs(i);
                    end
                    longspec(1000*(i-1)+1:1000*i) = Dnew; % Longspec might not produce a sensible spectrum, it is used to store spectra only
                end
                pydict = pyrunfile("python/predict_window3.py", "out_d", longwin=longspec); % Only runs once
                parastruct = struct(pydict); % Params are behind "s<i>p<k>", 'i' = no. split, 'k' = no. peak predicted
                for i=1:length(splits)-1
                    errs = zeros(1, numnets+1);
                    params = cell(1, numnets);
                    xerr = linspace(0,10,1000);
                    Dnew = longspec(1000*(i-1)+1:1000*i); % Might be scaled additionally
                    for k = 1:numnets
                        name = 's' + string(i) + 'p' + string(k);
                        pas = parastruct.(name); % Dynamic field name, no eval() needed
                        params{k} = double(pas); % Centers, heights, hws, lambdas
                        Derr = zeros(size(xerr));
                        for np= 1:k
                            Derr = Derr + gl(xerr, 1/2*params{k}(2*k+np), 10*params{k}(np), params{k}(k+np), 1/5*params{k}(3*k+np));
                        end
                        errs(k) = norm(Dnew-Derr);
                    end
                    % Find out how the extrapolated prediction compares
                    pind = find(tempMod.init.center(:,tempind) <= splits(i+1) & tempMod.init.center(:, tempind) >= splits(i)); % Finds indices of peaks in current range
                    Dex = zeros(size(Derr));
                    for np = 1:length(pind)
                        camt = tempMod.init.amt(pind(np),tempind);
                        ccent = 10*(tempMod.init.center(pind(np), tempind)-splits(i))/(splits(i+1)-splits(i));
                        cdist = max(-1, 10*tempMod.init.dist(pind(np), tempind)/(splits(i+1)-splits(i)));
                        chw  = 10*tempMod.init.hw(pind(np), tempind)/(splits(i+1)-splits(i));
                        cheight = tempMod.init.height(pind(np), tempind)*heightfacs(i);
                        clamb = tempMod.init.lambda(pind(np), tempind);
                        Dex = Dex + glmult(xerr, camt, ccent, cdist, chw, cheight, clamb, obj.BC);
                    end
                    exprederr = norm(Dnew-Dex);
                    errs(end) = min(norm(Dnew), exprederr);
                    if pps(i) == 0 || pps(i) > numnets
                        errs(end) = errs(end)*opt.ML.peakadv;
                    else
                        errs(pps(i)) = errs(pps(i))*opt.ML.peakadv; % Else let modif. best fit decide
                    end
                    [~, numpeak] = min(errs);
                    if numpeak ~= length(errs) % Then peaks are present, else they aren't, or every prediction is worse than not predicting.
                        left = splits(i);
                        right = splits(i+1);
                        sz = right-left;
                        temppara = params{numpeak};
                        temppara(1:numpeak) = temppara(1:numpeak)*sz + left; % Scale centers
                        temppara(numpeak+1:2*numpeak) = 1/heightfacs(i)*temppara(numpeak+1:2*numpeak); % Scale heights back
                        temppara(2*numpeak+1:3*numpeak) = sz/20*temppara(2*numpeak+1:3*numpeak); % Scale hws (1/2*1/10*sz)
                        temppara(3*numpeak+1:end) = 1/5*temppara(3*numpeak+1:end); % Scale lambdas
                        paracell{i} = temppara;
                    end
                end

                % After paracell has been obtained need to sort obtained parameters and peaks into next layer of tempmodel
                peakppm = zeros(1, sum(tempMod.init.amt(:,tempind))); % Assumed ppm-values of multiplet's subpeaks
                multipidx = peakppm; % Memory which peak belonged to which multiplet
                mnum = peakppm; % Memory where in multiplet current peak stands
                pind = 1;
                for p = 1:size(tempMod.init.amt, 1)
                    if tempMod.init.amt(p, tempind) == 1
                        peakppm(pind) = tempMod.init.center(p, tempind);
                        multipidx(pind) = p;
                        mnum(pind) = 1;
                        pind = pind + 1;
                    else
                        for pp = 1:tempMod.init.amt(p, tempind)
                            peakppm(pind) = tempMod.init.center(p, tempind) + (pp-1)*tempMod.init.dist(p, tempind);
                            multipidx(pind) = p;
                            mnum(pind) = pp;
                            pind = pind + 1;
                        end
                    end
                end
                tempstruct = tempMod.init; % Allowed to overwrite things in tempstruct
                fns = fieldnames(tempstruct);
                for n = [2,4,5,6] % Only need to overwrite certain fields in the struct
                    tempstruct.(fns{n}) = zeros(size(tempstruct.(fns{n})));
                end
                tempstruct.prednum = zeros(size(tempstruct.center,1),1);
                for ind = 1:length(paracell)
                    if ~isempty(paracell{ind})
                        logarr = min(peakppm>=splits(ind),peakppm<=splits(ind+1)); % Minimum represents logical AND
                        splitc = peakppm(logarr); 
                        splitinds = multipidx(logarr);
                        smnum = mnum(logarr);
                        paramat = reshape(paracell{ind}, [], 4)'; % [centers; heights; hws; lambdas]
                        [~, order] = sort(paramat(1,:));
                        paramat(:, :) = paramat(:, order); % Sorts matrix by first row (should be sorted, but with ML you never know)
                        im = length(splitc);
                        jm = size(paramat,2);
                        distmat = zeros(im, jm); % We calculate the distances of all predicted peaks to the extrapolated peaks and sort the nearest first
                        for iim = 1:im
                            for jjm = 1:jm
                                distmat(iim, jjm) = abs(splitc(iim)-paramat(1,jjm));
                            end
                        end
                        while im > 0 && jm > 0
                            [Ai, Aj] = find(distmat == min(distmat(:)),1);
                            distmat(Ai, :) = ones(size(distmat(Ai, :)))*inf;
                            distmat(:, Aj) = ones(size(distmat(:, Aj)))*inf;
                            im = im-1;
                            jm = jm-1;
                            nump = smnum(Ai);
                            Ai = splitinds(Ai);
                            % We take the average over all predictions of the corresponding multiplet parameters
                            if nump > 1
                                tempstruct.center(Ai, tempind) = tempstruct.center(Ai, tempind) + paramat(1, Aj) - tempstruct.dist(Ai, tempind)*(nump-1); % adjust centers to the fist peak of the multiplet
                                tempstruct.height(Ai, tempind) = tempstruct.height(Ai, tempind) + paramat(2, Aj)/obj.BC(tempstruct.amt(Ai, tempind), nump); % adjust heights
                            else
                                tempstruct.center(Ai, tempind) = tempstruct.center(Ai, tempind) + paramat(1, Aj);
                                tempstruct.height(Ai, tempind) = tempstruct.height(Ai, tempind) + paramat(2, Aj);
                            end
                            tempstruct.hw(Ai, tempind) = tempstruct.hw(Ai, tempind) + paramat(3, Aj);
                            tempstruct.lambda(Ai, tempind) = tempstruct.lambda(Ai, tempind) + paramat(4, Aj);
                            tempstruct.prednum(Ai) = tempstruct.prednum(Ai)+1;
                        end
                    end
                end
                % Prevents division by 0
                lvec = tempstruct.prednum > 0; % No changes in lvec, no need to change the init
                for n = [2 4 5 6] % Averaging and initializing tempMod
                    tempMod.init.(fns{n})(lvec, tempind) = tempstruct.(fns{n})(lvec, tempind)./tempstruct.prednum(lvec);
                end
                % Even further plots
                if opt.info == 2.2
                    disp(['Plot initial prediction of layer # ' num2str(Tsidx(j))])
                    d = D{Tsidx(j)};
                    x = X{Tsidx(j)};
                    MLmodel = zeros(size(x));
                    for p = 1:size(tempMod.init.center, 1)
                        MLmodel = MLmodel + glmult(x, tempMod.init.amt(p, tempind), tempMod.init.center(p, tempind), tempMod.init.dist(p, tempind), tempMod.init.hw(p, tempind), tempMod.init.height(p, tempind), tempMod.init.lambda(p, tempind), obj.BC);
                    end
                    figure(7);
                    set(gcf, 'Position', get(0, 'Screensize'));
                    if opt.epsinits
                        clf
                        hold on
                        plot(x,d,'k', 'LineWidth', 2)
                        plot(x,MLmodel, 'b', 'LineWidth', 2)
                    else
                        plot(x, d, 'k', x, MLmodel, 'b');
                    end
                    set(gca, 'XDir', 'reverse')
                    set(gcf, 'renderer', 'painters');
                    title('t1')
                    xlabel('x1')
                    xlim([xmin, xmax])
                    ylabel('y1')
                    set(gca, 'Fontsize', 24)
                    if opt.epsinits
                        ymax = max(max(d),max(MLmodel));
                        ylim([0,ymax])
                        set(gca, 'Fontsize', 48)
                        path = ['Bilder/Ints/abc128', num2str(Tsidx(j))];
                        pause(2)
                        print('-depsc', '-loose', path);
                        hold off
                    end
                    wm = wasserstein(x, MLmodel, d)*100/(x(end)-x(1)); % Normalized Wasserstein metric for one spectrum
                    frobm = norm(MLmodel - d, 'fro')*100/length(x); % MSE
                    disp(['Wasserstein-Metric of: ' num2str(wm) ' %'])
                    disp(['MSE of:' num2str(frobm) ' %'])
                end
                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                SpStatus(j) = 2;
            end
            if obj.options.edgeopt && opt.pred_layer_count < length(Tsidx) % Optimize first and last few layers extra
                for j=2:opt.pred_layer_count
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(1:j)])

                    optIDX = 1:j;

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getBoundedModel(optIDX, multidx);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                end
                for j=length(Tsidx)-1:-1:length(Tsidx)-opt.pred_layer_count+1
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(j:length(Tsidx))])

                    optIDX = j:length(Tsidx);

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getBoundedModel(optIDX, multidx);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                end
            end
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variant 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Initializing thin optimization process by setting optimization of centers and dists to false
        function obj = initElseOpt(obj)
            obj.model.act.center = false(size(obj.model.act.center));
            obj.model.act.dist = false(size(obj.model.act.dist));
            obj.model.init.center = obj.model.opt.center;
            obj.model.init.dist = obj.model.opt.dist;
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variant 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Finding vertical lines, at which no peaks crosses through, to split the optimization up
        function [splits, intidx] = findOptbounds(obj)
            X = obj.data.X;
            lb = zeros(1,obj.model.lenPeak);
            ub = lb;
            intidx = lb; % 'multidx(i)' = interval in which multiplet 'i' lives
            xmin = X{1}(1);
            xmax = X{1}(end);
            splits = [];

            % Finding unsplittable Intervals
            centers = obj.model.opt.center;
            amts = obj.model.opt.amt(:,1);
            dists = obj.model.opt.dist;
            for i=1:obj.model.lenPeak
                lb(i) = min(centers(i,:));
                ub(i) = max(centers(i,:)+(amts(i)-1)*dists(i,:));
            end
            [lb, order] = sort(lb, 'asc'); % Does the same as putting 1:k into 'intidx' and then ordering
            ub = ub(order);
            for i=1:length(intidx)
                intidx(i) = find(order == i);
            end
            lb(1) = xmin; % We don't want to split off to the left or right of all multiplets. This should be done preemptively
            ub(ub==max(ub)) = xmax;

            % Merging intervals
            splitdist = obj.options.L1.splitdist;
            ind = 1;
            while ub(ind) ~= xmax
                if lb(ind+1)-ub(ind)-splitdist < 0
                    % Merge intervals
                    ub(ind) = max(ub(ind), ub(ind+1)); % 'lb' already correct
                    lb(ind+1) = [];
                    ub(ind+1) = [];
                    intidx(intidx > ind) = intidx(intidx> ind)-1; % All numbers above 'ind' decrease by one
                else
                    % No merging needed, split found
                    splits(ind) = (ub(ind)+lb(ind+1))/2; %#ok<AGROW>  % Split in the middle
                    ind = ind+1;
                end
            end
            intidx(intidx > ind) = ind;
            splits = [xmin splits xmax];
        end
        %% Optimizing all variables occurring in certain bounds
        function obj = runBoundedOpt(obj, TsStartIdx, lb, ub, multidx)
            % 'TsStartIdx' is the index of 'Tsidx' to start with
            opt = obj.options;
            X = obj.data.X;
            D = obj.data.D;
            xmin = lb;
            xmax = ub;
            Tsidx = obj.model.Tsidx;
            ma = obj.data.ma;
            changed = {'center', 'dist'};

            % Limiting 'X' and 'D'
            lbind = find(X{1}>= lb, 1, 'first');
            ubind = find(X{1}<= ub, 1, 'last');

            for t = 1:length(X)
                X{t} = X{t}(lbind:ubind); % Not permanent
                D{t} = D{t}(lbind:ubind);
                ma(t) = max(abs(D{t}));
            end

            % SpStatus: 0 ~ unoptimized, 1 ~ upcoming, 2 ~ already optimized
            SpStatus = zeros(size(Tsidx));
            SpStatus(TsStartIdx) = 1;

            % Running initial optimization
            disp(['Run initial: Selected layer # ' num2str(TsStartIdx)])
            optIDX = TsStartIdx;
            tempMod = obj.model.getBoundedModel(optIDX, multidx); % Needs creating
            tempMod.init = tempMod.opt; % Needs corrections
            tempMod.lower.center = max(tempMod.lower.center, lb*ones(size(tempMod.lower.center))); % If somehow stricter limits exist, keep those
            tempMod.upper.center = min(tempMod.upper.center, ub*ones(size(tempMod.upper.center)));

            opttime = tempMod.optimizeModel(opt, X(Tsidx(optIDX)), D(Tsidx(optIDX)), ma(Tsidx(optIDX)), obj.spath, obj.BC); % 'X(inds)' etc. gives cell array, not its contents
            disp(['Computation time: ' num2str(opttime) 's'])
            % Transferring results back to original model
            obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx); % Needs adapting
            SpStatus(optIDX) = 2;

            % Running forward
            for j = TsStartIdx +1 : length(Tsidx)
                disp(['Run forward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = max(TsStartIdx,j-opt.pred_layer_count);
                ide = j;
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getBoundedModel(optIDX, multidx);
                tempMod.init = tempMod.opt;
                tempMod.lower.center = lb*ones(size(tempMod.lower.center));
                tempMod.upper.center = ub*ones(size(tempMod.upper.center));
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax, changed);  % Harmless

                % Look, plots
                if opt.info == 2.2
                    tempind = find(SpStatus(optIDX)==1);
                    disp(['Plot initial prediction of layer # ' num2str(Tsidx(j))])
                    d = D{Tsidx(j)};
                    x = X{Tsidx(j)};
                    MLmodel = zeros(size(x));
                    for p = 1:size(tempMod.init.center, 1)
                        MLmodel = MLmodel + glmult(x, tempMod.init.amt(p, tempind), tempMod.init.center(p, tempind), tempMod.init.dist(p, tempind), tempMod.init.hw(p, tempind), tempMod.init.height(p, tempind), tempMod.init.lambda(p, tempind), obj.BC);
                    end
                    figure(7);
                    set(gcf, 'Position', get(0, 'Screensize'));
                    if opt.epsinits
                        clf
                        hold on
                        plot(x,d,'k', 'LineWidth', 2)
                        plot(x,MLmodel, 'b', 'LineWidth', 2)
                    else
                        plot(x, d, 'k', x, MLmodel, 'b');
                    end
                    set(gca, 'XDir', 'reverse')
                    set(gcf, 'renderer', 'painters');
                    title('t1')
                    xlabel('x1')
                    xlim([xmin, xmax])
                    ylabel('y1')
                    set(gca, 'Fontsize', 24)
                    if opt.epsinits
                        ymax = max(max(d),max(MLmodel));
                        ylim([0,ymax])
                        set(gca, 'Fontsize', 48)
                        path = ['Bilder/Ints/abc129', num2str(Tsidx(j))];
                        pause(2)
                        print('-depsc', '-loose', path);
                        hold off
                    end
                    wm = wasserstein(x, MLmodel, d)*100/(x(end)-x(1)); % Normalized Wasserstein metric for one spectrum
                    frobm = norm(MLmodel - d, 'fro')*100/length(x); % MSE
                    disp(['Wasserstein-Metric of: ' num2str(wm) ' %'])
                    disp(['MSE of:' num2str(frobm) ' %'])
                end
                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                SpStatus(j) = 2;
            end

            % Running backward
            for j = TsStartIdx -1:-1:1
                disp(['Run backward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = j;
                ide = min(length(SpStatus),j+opt.pred_layer_count); 
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getBoundedModel(optIDX, multidx);
                tempMod.init = tempMod.opt;
                tempMod.lower.center = lb*ones(size(tempMod.lower.center));
                tempMod.upper.center = ub*ones(size(tempMod.upper.center));
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax, changed);

                % Look, more plots
                if opt.info == 2.2
                    tempind = find(SpStatus(optIDX)==1);
                    disp(['Plot initial prediction of layer # ' num2str(Tsidx(j))])
                    d = D{Tsidx(j)};
                    x = X{Tsidx(j)};
                    MLmodel = zeros(size(x));
                    for p = 1:size(tempMod.init.center, 1)
                        MLmodel = MLmodel + glmult(x, tempMod.init.amt(p, tempind), tempMod.init.center(p, tempind), tempMod.init.dist(p, tempind), tempMod.init.hw(p, tempind), tempMod.init.height(p, tempind), tempMod.init.lambda(p, tempind), obj.BC);
                    end
                    figure(7);
                    set(gcf, 'Position', get(0, 'Screensize'));
                    if opt.epsinits
                        clf
                        hold on
                        plot(x,d,'k', 'LineWidth', 2)
                        plot(x,MLmodel, 'b', 'LineWidth', 2)
                    else
                        plot(x, d, 'k', x, MLmodel, 'b');
                    end
                    set(gca, 'XDir', 'reverse')
                    set(gcf, 'renderer', 'painters');
                    title('t1')
                    xlabel('x1')
                    xlim([xmin, xmax])
                    ylabel('y1')
                    set(gca, 'Fontsize', 24)
                    if opt.epsinits
                        ymax = max(max(d),max(MLmodel));
                        ylim([0,ymax])
                        set(gca, 'Fontsize', 48)
                        path = ['Bilder/Ints/abc130', num2str(Tsidx(j))];
                        pause(2)
                        print('-depsc', '-loose', path);
                        hold off
                    end
                    wm = wasserstein(x, MLmodel, d)*100/(x(end)-x(1)); % Normalized Wasserstein metric for one spectrum
                    frobm = norm(MLmodel - d, 'fro')*100/length(x); % MSE
                    disp(['Wasserstein-Metric of: ' num2str(wm) ' %'])
                    disp(['MSE of:' num2str(frobm) ' %'])
                end
                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                SpStatus(j) = 2;
            end

            if obj.options.edgeopt && opt.pred_layer_count < length(Tsidx) % Optimize first and last few layers extra
                for j=2:opt.pred_layer_count
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(1:j)])

                    optIDX = 1:j;

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getBoundedModel(optIDX, multidx);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                end
                for j=length(Tsidx)-1:-1:length(Tsidx)-opt.pred_layer_count+1
                    disp(['Run edge-time-layers: Selected layer(s) # ' num2str(j:length(Tsidx))])

                    optIDX = j:length(Tsidx);

                    % Constructing the corresponding reduced model
                    tempMod = obj.model.getBoundedModel(optIDX, multidx);
                    tempMod.init = tempMod.opt;

                    % Actual optimization process
                    range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                    opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                    disp(['Computation time: ' num2str(opttime) 's'])
                    % Transferring results back to original model
                    obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                end
            end
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variant 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Optimizing all variables occurring in certain bounds
        function obj = runSausageOpt(obj, TsStartIdx, multidx)
            % 'TsStartIdx' is the index of 'Tsidx' to start with
            opt = obj.options;
            X = obj.data.X;
            D = obj.data.D;
            Tsidx = obj.model.Tsidx;
            T = obj.data.T;
            ma = obj.data.ma;
            changed = {'center', 'dist'};
            centers = obj.model.opt.center(multidx, :);
            amt = obj.model.opt.amt(multidx, 1);
            dists = obj.model.opt.dist(multidx, :);

            centerfull = pchip(T(Tsidx), centers, T);
            distfull = pchip(T(Tsidx), dists, T);
            SausageShift = opt.L1opt.sshift;

            xmin = X{1}(1);
            xmax = X{1}(end);

            % Limiting 'X' and 'D'
            lb = max(centerfull-SausageShift, xmin);
            ub = min(centerfull+(amt-1)*distfull+SausageShift, xmax);

            for t = 1:length(X)
                lbind = find(X{t}>= lb(t), 1, 'first');
                ubind = find(X{t}<= ub(t), 1, 'last');
                X{t} = X{t}(lbind:ubind); % Nothing is permanent
                D{t} = D{t}(lbind:ubind);
                ma(t) = max(abs(D{t}));
            end

            % SpStatus: 0 ~ unoptimized, 1 ~ upcoming, 2 ~ already optimized
            SpStatus = zeros(size(Tsidx));
            SpStatus(TsStartIdx) = 1;

            % Running initial optimization
            disp(['Run initial: Selected layer # ' num2str(TsStartIdx)])
            optIDX = TsStartIdx;
            tempMod = obj.model.getBoundedModel(optIDX, multidx); % Same function should work
            tempMod.init = tempMod.opt; % Needs corrections
            tempMod.lower.center = max(tempMod.opt.center - SausageShift, xmin);
            tempMod.upper.center = min(tempMod.opt.center + (amt-1)*tempMod.opt.dist + SausageShift, xmax);

            opttime = tempMod.optimizeModel(opt, X(Tsidx(optIDX)), D(Tsidx(optIDX)), ma(Tsidx(optIDX)), obj.spath, obj.BC); % X(inds) etc. gives cell array, not its contents
            disp(['Computation time: ' num2str(opttime) 's'])
            % Transferring results back to original model
            obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx); % Has been adapted
            SpStatus(optIDX) = 2;

            % Running forward
            for j = TsStartIdx +1 : length(Tsidx)
                disp(['Run forward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = max(TsStartIdx,j-opt.pred_layer_count);
                ide = j;
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getBoundedModel(optIDX, multidx);
                tempMod.init = tempMod.opt;
                tempMod.lower.center = max(tempMod.opt.center - SausageShift, xmin);
                tempMod.upper.center = min(tempMod.opt.center + (amt-1)*tempMod.opt.dist + SausageShift, xmax);

                tempMod.interpPara(SpStatus(optIDX), xmin, xmax, changed);

                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                SpStatus(j) = 2;
            end

            % Running backward
            for j = TsStartIdx -1:-1:1
                disp(['Run backward in time: Selected layer # ' num2str(j)])
                SpStatus(j) = 1;
                filt = zeros(size(SpStatus));
                ids = j;
                ide = min(length(SpStatus),j+opt.pred_layer_count);
                filt(ids:ide) = 1;
                optIDX = find(filt);

                % Constructing the corresponding reduced model
                tempMod = obj.model.getBoundedModel(optIDX, multidx);
                tempMod.init = tempMod.opt;
                tempMod.lower.center = max(tempMod.opt.center - SausageShift, xmin);
                tempMod.upper.center = min(tempMod.opt.center + (amt-1)*tempMod.opt.dist + SausageShift, xmax);
                tempMod.interpPara(SpStatus(optIDX), xmin, xmax, changed);

                % Actual optimization process
                range = min(Tsidx(optIDX)):max(Tsidx(optIDX));
                opttime = tempMod.optimizeModel(opt, X(range), D(range), ma(range), obj.spath, obj.BC);
                disp(['Computation time: ' num2str(opttime) 's'])
                % Transferring results back to original model
                obj.model.insertBoundedPara(tempMod.opt, optIDX, multidx);
                SpStatus(j) = 2;
            end
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variant 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Calling single-layer-analysis function of current model
        function obj = analyzeSingLayer(obj, specidx)
            Xinit = obj.data.X{specidx}(:);
            Dinit = obj.data.D{specidx}(:);
            obj.model.analyzeSingLayer(Xinit, Dinit, specidx, obj.options, obj.BC);
        end
    end
end