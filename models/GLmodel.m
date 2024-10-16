classdef GLmodel < model
    %  Convex combination of Gaussian and Lorentzian models, inherits from class model. Only multiplets.

    properties
        mname = 'Gauss-Lorentzmodel'% Predefined name of model - don't change
        lenPeak                     % Number of Multiplets present
        lenPara                     % #Parameters * #SplineLayers * #Multiplets
        Tfull                       % Full time vector
        Ts                          % Spline knot time vector
        Tsidx                       % Indices of Ts (Tfull(Tsidx) = Ts)

        init                        % Struct init.amt/center/dist/hw/height/lambda
        % .amt = amount of peaks in multiplet           - matrix lenPeak x length(Ts)
        % .center = center values                       - matrix lenPeak x length(Ts)
        % .dist  = distance between peaks of multiplet  - matrix lenPeak x length(Ts)
        % .hw  = half width values                      - matrix lenPeak x length(Ts)
        % .height = height values                       - matrix lenPeak x length(Ts)
        % .lambda = convex combination factor           - matrix lenPeak x length(Ts)
        lower                       % Lower opt. bounds - struct
        upper                       % Upper opt. bounds - struct
        act                         % 'Active' Boolean  - struct
        opt                         % Optimized values  - struct

        scaleM                      % Scaling matrix for parameters
        scaleMi                     % Inverted scaling matrix

        ParaName = {'amt', 'center', 'dist', 'hw', 'height', 'lambda'} % Don't change!

    end

    methods
        %% Constructor
        function obj = GLmodel(Tfull, init, lower, upper, act, Tsidx)
            obj.Tfull = Tfull;

            if nargin == 1     % Empty models will be important, see next method
                obj.Ts = [];

                obj.init = [];
                obj.lower = [];
                obj.upper = [];
                obj.opt = [];
                obj.act = [];

                obj.lenPeak = 0;
                obj.lenPara = 0;
            else
                obj.Tsidx = Tsidx;
                obj.Ts = obj.Tfull(Tsidx);

                obj.init = init;
                obj.lower = lower;
                obj.upper = upper;
                obj.opt = init ;
                obj.act = act;

                obj.lenPeak = size(init.center, 1);
                obj.lenPara = length(obj.ParaName)*obj.lenPeak*length(obj.Ts);
            end
        end

        %% Getting reduced model of this existing model (of a few timelayers of the SplineLayers)
        function redmod = getReducedModel(obj, redTsidxIdx)
            redTsidx = obj.Tsidx(redTsidxIdx);

            redmod = GLmodel(obj.Tfull(min(redTsidx):max(redTsidx)));

            redmod.Ts = obj.Tfull(redTsidx);
            redmod.Tsidx = redTsidx - min(redTsidx)+1; % In the reduced model the timelayers start at Tfull(min(redTsidx))
            lenTs = length(redmod.Ts);
            redmod.lenPeak = obj.lenPeak;
            redmod.lenPara = length(redmod.ParaName)*redmod.lenPeak*lenTs;

            for k=1:length(redmod.ParaName)   % Eval prevents writing 6 blocks of very similar code.
                eval(['redmod.init.' obj.ParaName{k} '=obj.init.' obj.ParaName{k} '(:, redTsidxIdx);'])
                eval(['redmod.lower.' obj.ParaName{k} '=obj.lower.' obj.ParaName{k} '(:, redTsidxIdx);'])
                eval(['redmod.upper.' obj.ParaName{k} '=obj.upper.' obj.ParaName{k} '(:, redTsidxIdx);'])
                eval(['redmod.act.' obj.ParaName{k} '=obj.act.' obj.ParaName{k} '(:, redTsidxIdx);'])
                eval(['redmod.opt.' obj.ParaName{k} '=obj.opt.' obj.ParaName{k} '(:, redTsidxIdx);'])
            end

            for i=1:length(obj.constraints)
                redmod.addConstraint(obj.constraints{i}.type, obj.constraints{i}.weight, obj.constraints{i}.var, obj.constraints{i}.Pidx);
            end
        end

        %% Inserting parameters in structform into certain special timelayers
        function insertPara(obj,inp_opt,redTsidxIdx)  %#ok<INUSD>
            for k=1:length(obj.ParaName)
                eval(['obj.opt.' obj.ParaName{k} '(:,redTsidxIdx)=inp_opt.' obj.ParaName{k} ';'])
            end
        end

        %% Extrapolating parameters to next SplineLayers (also needs some interpolation, i guess)
        % Status represents the opimization status in the optimization step of all time layers
        % 0 ~ unoptimized, 1 ~ needs optimization now, 2 ~ has already been optimized, is used to extrapolate
        function interpPara(obj, status, xmin, xmax, changed)
            idxint = find(status == 1);
            tint = obj.Ts(idxint); %#ok<NASGU>
            idxbase = find(status == 2);
            tbase = obj.Ts(idxbase); %#ok<NASGU>
            actflag = 1;
            changeflag = 0;

            if nargin <= 4
                changed = {'foobar'};
            end

            for i=1:obj.lenPeak
                for k = 2:length(obj.ParaName) % No need to change 'amt'
                    eval(['actflag = obj.act.' obj.ParaName{k} '(i,idxint);']) % Do we need to interpolate (only relevant for split optimization)
                    if max(strcmp(changed, obj.ParaName{k}))
                        changeflag = 1;
                    end
                    interpflag = actflag && (~changeflag);
                    if interpflag
                        eval(['ybase = obj.init.' obj.ParaName{k} '(i,idxbase);']) % 'ybase' only in use for one parameter at a time
                        if length(ybase) == 1
                            eval(['obj.init.' obj.ParaName{k} '(i,idxint) = obj.init.' obj.ParaName{k} '(i,idxbase);'])
                        else
                            eval(strjoin(['obj.init.' obj.ParaName{k} "(i,idxint) = interp1(tbase,ybase ,tint,'linear','extrap');"],''))
                        end
                    end
                end
                obj.init.amt(i, idxint) = obj.init.amt(i, idxbase(1)); % 'amt' stays constant everywhere
                if obj.init.amt(i,idxint) == 1
                    obj.init.dist(i,idxint) = -1;
                end
                obj.init.height(i, idxint) = abs(obj.init.height(i, idxint));
                obj.init.hw(i, idxint) = abs(obj.init.hw(i, idxint));
                obj.init.lambda(i, idxint) = min(max(obj.init.lambda(i,idxint),0),1); % Always start between 0, 1
                obj.init.center(i, idxint) = min(xmax - 0.05, max(xmin + 0.05, obj.init.center(i, idxint))); % Being at the edge is disadvantageous
            end
        end

        %% Initializing all parameters on SplineLayers
        function obj = initializeConst(obj, initLayer, Tsidx, options, ma, X)
            %'initLayer' contains information and values for first time layer in a struct
            obj.Tsidx = Tsidx;
            obj.Ts = obj.Tfull(Tsidx);
            lenTs = length(obj.Ts);
            obj.lenPeak = length(initLayer.center);
            obj.lenPara = length(obj.ParaName)*obj.lenPeak*lenTs;

            for k=2:length(obj.ParaName) % 'amt' behaves differently
                eval(['obj.init.' obj.ParaName{k} ' = initLayer.' obj.ParaName{k} '(:)*ones(1,lenTs);']) % Produces a matrix
                eval(['obj.act.' obj.ParaName{k} ' =  true(obj.lenPeak, lenTs);'])
            end

            obj.init.amt = initLayer.amt(:)*ones(1,lenTs);
            obj.act.amt  = false(obj.lenPeak, lenTs);
            for i=1:obj.lenPeak
                if obj.init.amt(i,1) == 1
                    obj.init.dist(i,:) = -1;
                    obj.act.dist(i,:) = false; % Yes, that works. Don't optimize distance if there's none to optimize.
                end
                if ~options.L1opt.optdist
                    obj.act.dist(i,:) = false; % Maybe it is not desirable to optimize distance
                end
            end
            obj.lower.center = zeros(size(obj.init.center)); % Otherwise dimension goes awry in following loop
            obj.upper.center = zeros(size(obj.init.center));
            obj.upper.height = zeros(size(obj.init.height));

            for i=1:lenTs
                obj.lower.center(:,i) = min(X{Tsidx(i)});   % Limits centers
                obj.upper.center(:,i) = max(X{Tsidx(i)});
                obj.upper.height(:,i) = 1.5*ma(Tsidx(i));   % Peaks shouldn't be much higher than the data, no need to optimize further
            end
            obj.lower.amt = initLayer.amt(:)*ones(1,lenTs); % Never should change, implement same upper and lower bounds
            obj.lower.dist = zeros(obj.lenPeak, lenTs); % If 'act' is 0 for one 'dist', bounds don't matter
            obj.lower.hw = zeros(obj.lenPeak, lenTs);
            obj.lower.height = zeros(obj.lenPeak, lenTs);
            obj.lower.lambda = options.L1opt.minlamb*ones(obj.lenPeak, lenTs);

            obj.upper.amt = initLayer.amt(:)*ones(1,lenTs); % Same reason
            obj.upper.dist = options.L1opt.maxdist*ones(obj.lenPeak, lenTs);
            obj.upper.hw = options.L1opt.maxhw*ones(obj.lenPeak, lenTs);
            obj.upper.lambda = options.L1opt.maxlamb*ones(obj.lenPeak, lenTs);

            obj.opt = obj.init;
        end

        %% Optimizing the model's first layer
        function initLayer = analyzeFirstLayer(obj, params, Xinit, Dinit, options, spath, BC)
            % After finding all 'centeridx' & 'amt' values and suggested
            % 'distidx' (params), optimize 'dist', 'hw', 'height' and 
            % 'centeridx' in ascending order, rest corresponding to 
            % 'centeridx' respectively
            L1opt = options.L1opt;
            bw = options.bwdth;
            adapt = options.L1opt.adapt;

            obj.lenPeak         = length(params.centeridx);
            initLayer.center    = Xinit(params.centeridx);
            initLayer.center    = initLayer.center(:); % Produces a column vector
            initLayer.amt       = params.amt(:);
            initLayer.dist      = zeros(size(initLayer.center));
            initLayer.hw        = zeros(size(initLayer.center));
            initLayer.height    = zeros(size(initLayer.center));
            initLayer.lambda    = zeros(size(initLayer.center));

            lefts = zeros(length(params.centeridx),1); % Stores indices of local optimization boundaries
            rights = zeros(length(params.centeridx),1);

            for j = 1:obj.lenPeak % Multiplets are easy
                if initLayer.amt(j) ~= 1
                    lefts(j) = max(1,params.centeridx(j) - max(bw, params.distidx(j))); % To the left of the whole multiplet
                    rights(j) = min(length(Xinit),params.centeridx(j) + max(bw*params.amt(j), params.distidx(j)*params.amt(j))); % To the right
                end
            end

            if (obj.lenPeak == 1 && initLayer.amt(1) == 1) % A singular peak is easy
                lefts(1) = max(1,params.centeridx(j) - adapt);
                rights(1) = min(length(Xinit),params.centeridx(j) + adapt);
            else
                for j=find(lefts==0)' % Remaining singlets:
                    if j == 1
                        lefts(j) = max(1, params.centeridx(j)-adapt); % Don't overshoot, and on right side don't overshoot 'adapt' or next peak as well
                        rights(j)= min([length(Xinit), params.centeridx(j)+adapt, ceil((params.centeridx(j)+params.centeridx(j+1))/2)]);
                    elseif j == obj.lenPeak
                        if rights(j-1) < params.centeridx(j) % Don't overshoot 'rights(j-1)' as well, since it denotes the end of a multiplet
                            lefts(j) = max([1, params.centeridx(j)-adapt, rights(j-1)]);
                            rights(j) = min(length(Xinit), params.centeridx(j)+adapt);
                        else  % Is only the case if the current peak is inside of another multiplet
                            lefts(j) = max(1, params.centeridx(j)- bw); % We should have 'bw' << 'adapt'
                            rights(j) = min(length(Xinit), params.centeridx(j)+bw);
                        end

                    else
                        if rights(j-1) < params.centeridx(j)
                            lefts(j) = max([1, params.centeridx(j)-adapt, rights(j-1)]);
                            rights(j)= min([length(Xinit), params.centeridx(j)+adapt, ceil((params.centeridx(j)+params.centeridx(j+1))/2)]);
                        else
                            lefts(j) = max(1, params.centeridx(j)- bw);
                            rights(j) = min(length(Xinit), params.centeridx(j)+bw);
                        end

                    end
                end
            end


            optopt = optimset('Display','iter','MaxIter',options.Init.MaxIter,'MaxFunEvals',options.Init.MaxFunEvals,'TolFun',options.Init.TolFun,'TolX',options.Init.TolX);

            % Optimize 'dist', 'hw', 'height' and 'lambda# at once, first for each peak individually
            distinit = zeros(obj.lenPeak, 1);
            ldistinit = distinit;
            udistinit = distinit;
            C = zeros(obj.lenPeak,1);
            for i=1:obj.lenPeak
                distinit(i) = abs(Xinit(params.centeridx(i))-Xinit(params.centeridx(i)+params.distidx(i))); % Set for all peaks, only used for true multiplets
                ldistinit(i) = abs(Xinit(params.centeridx(i))-Xinit(params.centeridx(i)+params.distidx(i)-options.L1.multtol));
                udistinit(i) = abs(Xinit(params.centeridx(i))-Xinit(params.centeridx(i)+params.distidx(i)+options.L1.multtol));
                C(i) = 0.5*Dinit(params.centeridx(i))*nchoosek(params.amt(i)-1, ceil((params.amt(i)-1)/2));
            end
            ainit = 0.01;
            lambinit = 0.1;
            LB_1 = [0,0,L1opt.minlamb];  % If amt == 1
            UB_1 = [L1opt.maxhw, inf, L1opt.maxlamb];
            LB_2 = [0,0,0, L1opt.minlamb]; % Else
            UB_2 = [L1opt.maxdist, L1opt.maxhw, inf,L1opt.maxlamb];

            for i = 1:obj.lenPeak
                lb = lefts(i);
                rb = rights(i);
                amt = params.amt(i);
                center = Xinit(params.centeridx(i));
                if amt == 1
                    smol_fun = @ (para) (Dinit(lb:rb) -glmult(Xinit(lb:rb), amt, center, -1, para(1), para(2), para(3), BC)); % errorfunction is difference from data to model
                    [para_opt, ~] = lsqnonlin(smol_fun, [ainit, C(i), lambinit], LB_1, UB_1, optopt);

                    initLayer.dist(i) = -1; % To indicate there's nothing to be optimized
                    initLayer.hw(i) = para_opt(1);
                    initLayer.height(i) = para_opt(2);
                    initLayer.lambda(i) = para_opt(3);
                else
                    LB_2(1) = ldistinit(i);
                    UB_2(1) = udistinit(i);
                    smol_fun = @ (para) (Dinit(lb:rb) -glmult(Xinit(lb:rb), amt, center, para(1), para(2), para(3), para(4), BC)); % errorfunction is difference from data to model
                    [para_opt, ~] = lsqnonlin(smol_fun, [distinit(i), ainit, C(i), lambinit], LB_2, UB_2, optopt);

                    initLayer.dist(i) = para_opt(1); % Otherwise there's the optimized value
                    initLayer.hw(i) = para_opt(2);
                    initLayer.height(i) = para_opt(3);
                    initLayer.lambda(i) = para_opt(4);
                end
            end
            % ~ Nice plots x5
            if options.info >=2 % Niche plots, again.
                figure;
                if options.info == 2.1
                    str = 'b';
                    hold on
                    col = [0.4660 0.6740 0.3];
                    lw = 1.5;
                else
                    str = 'k';
                    col = [0 0 1];
                    lw = 1;
                end
                plot(Xinit, Dinit, str);
                set(gca, 'XDir', 'reverse');
                hold on
                for i=1:obj.lenPeak
                    plot(Xinit, glmult(Xinit, initLayer.amt(i), initLayer.center(i), initLayer.dist(i), initLayer.hw(i), initLayer.height(i), initLayer.lambda(i), BC), 'r--')
                    plot(Xinit([lefts(i),lefts(i)]), [0, initLayer.height(i)], 'color', col, 'Linewidth', lw)
                    plot(Xinit([rights(i), rights(i)]), [0, initLayer.height(i)], 'color', col, 'Linewidth', lw)
                end
            end

            clower.amt       = initLayer.amt;
            cupper.amt       = initLayer.amt;
            cact.amt         = false(size(initLayer.amt)); %Don't optimize 'amt'.

            clower.center    = zeros(size(initLayer.center)); % Only temporary, no worries
            cupper.center    = zeros(size(initLayer.center));
            cact.center      = true(size(initLayer.center));

            clower.dist      = zeros(size(initLayer.dist));
            cupper.dist      = options.L1opt.maxdist*ones(size(initLayer.dist));
            cact.dist        = true(size(initLayer.dist)); % Only temporarily fully active

            clower.hw        = zeros(size(initLayer.hw));
            cupper.hw        = L1opt.maxhw*ones(size(initLayer.hw));
            cact.hw          = true(size(initLayer.hw));

            clower.height    = zeros(size(initLayer.height));
            cupper.height    = inf*ones(size(initLayer.height));
            cact.height      = true(size(initLayer.height));

            clower.lambda    = L1opt.minlamb*ones(size(initLayer.lambda));
            cupper.lambda    = L1opt.maxlamb*ones(size(initLayer.lambda));
            cact.lambda      = true(size(initLayer.lambda));

            for i=1:obj.lenPeak % See, i told you not to worry
                clower.center(i) = Xinit(max(params.centeridx(i)-L1opt.maxcenter,1));
                cupper.center(i) = Xinit(min(params.centeridx(i)+L1opt.maxcenter,length(Xinit)));
                if initLayer.amt(i) == 1
                    cact.dist(i,:) = false; % Again, this works
                end
            end

            model = GLmodel(obj.Tfull(1), initLayer, clower, cupper, cact, 1); % Creates auxiliary model of a single time layer to optimize
            model.optimizeModel(options, {Xinit}, {Dinit}, max(Dinit), spath, BC);

            initLayer = model.opt;

            % ~ Nice plots x6
            if options.info >=2
                figure;
                if options.info == 2.1
                    hold on
                end
                plot(Xinit, Dinit, str);
                set(gca, 'XDir', 'reverse');
                hold on
                for i=1:obj.lenPeak
                    plot(Xinit, glmult(Xinit, initLayer.amt(i), initLayer.center(i), initLayer.dist(i), initLayer.hw(i), initLayer.height(i), initLayer.lambda(i), BC), 'r--')
                end
            end
        end

        %% Calculating the model
        function Y = MOD(obj, para, X, BC) % Output in cell-format 'Y{i}(j)' has the datapoint 'j' in timelayer 'i'
            para = obj.scaleM*para; % Parametres are scaled 

            paraR = reshape(para,obj.lenPeak, length(obj.Ts), length(obj.ParaName)); % Shape matters, makes parameters more easily accessible

            amtfull     = zeros(obj.lenPeak, length(obj.Tfull));
            centerfull  = amtfull;
            distfull    = amtfull;
            hwfull      = amtfull;
            heightfull 	= amtfull;
            lambdafull  = amtfull;

            if length(obj.Ts) == 1 % Handles differently, doesn't need spline interpolation
                amtfull(:,1)    = round(paraR(:,1,1),0); % If 'scaleM' and 'scaleMi' destroy the integer
                centerfull(:,1) = paraR(:,1,2);
                distfull(:,1)   = paraR(:,1,3);
                hwfull(:,1)     = paraR(:,1,4);
                heightfull(:,1) = paraR(:,1,5);
                lambdafull(:,1) = paraR(:,1,6);

                Y = cell(1);
                Y{1} = zeros(size(X{1}));
                for i=1:obj.lenPeak
                    Y{1} = Y{1} + glmult(X{1}, amtfull(i,1), centerfull(i,1), distfull(i,1), hwfull(i,1), heightfull(i,1), lambdafull(i,1), BC);
                end
            else
                for i=1:obj.lenPeak
                    amtfull(i,:)    = paraR(i,1,1)*ones(size(amtfull(i,:))); % No risk, i don't trust pchip to always output integers
                    centerfull(i,:) = pchip(obj.Ts, paraR(i,:,2), obj.Tfull);
                    distfull(i,:)   = pchip(obj.Ts, paraR(i,:,3), obj.Tfull); % Doesn't need to strictly be an integer, if amt == 1 it's effectively unused
                    hwfull(i,:)     = pchip(obj.Ts, paraR(i,:,4), obj.Tfull);
                    heightfull(i,:) = pchip(obj.Ts, paraR(i,:,5), obj.Tfull);
                    lambdafull(i,:) = pchip(obj.Ts, paraR(i,:,6), obj.Tfull);
                end
                Y = cell(length(obj.Tfull),1);
                for j=1:length(obj.Tfull)
                    Y{j} = zeros(size(X{j}));
                    for i=1:obj.lenPeak
                        Y{j} = Y{j} + glmult(X{j}, amtfull(i,j), centerfull(i,j), distfull(i,j), hwfull(i,j), heightfull(i,j), lambdafull(i,j), BC);
                    end
                end
            end
        end

        %% Plotting the model (used after calculations)
        function obj = plot(obj, Spl_obj)
            figure;
            if nargin == 2
                data = Spl_obj.data;
                % Checks if common plot functions can be used
                b=1;
                for i=1:length(data.T)-1
                    if length(data.X{i})==length(data.X{i+1})
                        if max(abs(data.X{i}-data.X{i+1}))~=0
                            b=0;
                        end
                    else
                        b=0;
                    end
                end


                col1=[0    0.4470    0.7410];       % Colour scheme
                col2=[0.8500    0.3250    0.0980];

                if b
                    Dm = zeros(length(data.X{1}),length(data.T));
                    for j=1:length(data.T)
                        Dm(:,j) = data.D{j};
                    end

                    surf(data.T, data.X{1}, Dm); % Underlying plot below the red-ish peak movement lines
                    hold on
                    shading interp
                    view([-60 45]);
                else
                    for j=1:length(data.T)
                        l=j/length(data.T);
                        plot3(data.X{j}, data.T(j)*ones(size(data.X{j})), data.D{j}, 'color', l*col2 + (1-l)*col1);
                        hold on
                    end
                end
            end
            for i=1:obj.lenPeak % Plot the peaks of all multiplets: the further from the center peak, the lighter the shade of red
                for j=1:obj.opt.amt(i,1)
                    heightcoeff = nchoosek(obj.opt.amt(i,1)-1,j-1)/nchoosek(obj.opt.amt(i,1)-1,ceil((obj.opt.amt(i,1)-1)/2));
                    plot3(obj.Tfull, pchip(obj.Ts, obj.opt.center(i,:) + (j-1)*obj.opt.dist(i,:), obj.Tfull), pchip(obj.Ts, obj.opt.height(i,:), obj.Tfull)*heightcoeff, 'color', [1 0 0 (1-abs((obj.opt.amt(i,1)+1)/2 - j)/obj.opt.amt(i,1))], 'linewidth', 2);
                    hold on
                end
                plot3(obj.Ts, obj.opt.center(i,:) + ((obj.opt.amt(i,1)-1)/2)*obj.opt.dist(i,:), obj.opt.height(i,:), 'r*', 'linewidth', 2);
            end
            % Plot the initial values on all layers with 'bo's
            for i=1:obj.lenPeak
                plot3(obj.Ts, obj.init.center(i,:) + ((obj.init.amt(i,1)-1)/2)*obj.init.dist(i,:), obj.init.height(i,:), 'bo', 'linewidth',2);
            end
            axis tight
        end

        %% Optimizing the model (not just initially)
        function comptime = optimizeModel(obj, opt, X, D, ma, spath, BC)
            tic;
            parainit = obj.para2vec(obj.init); % Index magic
            paraLB = obj.para2vec(obj.lower);
            paraUB = obj.para2vec(obj.upper);
            paraAct = obj.para2vec(obj.act);

            nr_runs = length(opt.scheme);
            for i=1:nr_runs
                disp(['RUN ' num2str(i)]);

                parainitS = abs(parainit);
                parainitS(parainitS < 1E-6) = 1; % Prevents ill-conditioned scale-matrix
                obj.scaleM = diag(parainitS);
                obj.scaleMi = inv(obj.scaleM);
                parainit = obj.scaleMi*parainit; % Scaling might interfere with amt being an integer, amt is being rounded in obj.MOD
                paraLB = obj.scaleMi*paraLB;
                paraUB = obj.scaleMi*paraUB;

                % Give options to lsqnonlin
                if length(X) == 1
                    opt.TolX = opt.Init.TolX; % Structs are called by value, so changeing the struct here does not change the struct globally
                    opt.TolFun = opt.Init.TolFun;
                    opt.MaxFunEvals = opt.Init.MaxFunEvals;
                    opt.MaxIter = opt.Init.MaxIter;
                end
                switch opt.precon
                    case 1
                        optopt = optimoptions(@lsqnonlin, 'Display', 'Iter', 'MaxIter', opt.MaxIter, 'MaxFunEvals', opt.MaxFunEvals, 'TolFun', opt.TolFun, 'TolX', opt.TolX, 'plotfcn', @optimplotx);
                    case 2
                        optopt = optimoptions(@lsqnonlin, 'Display', 'Iter', 'MaxIter', opt.MaxIter, 'MaxFunEvals', opt.MaxFunEvals, 'TolFun', opt.TolFun, 'TolX', opt.TolX, 'plotfcn', @optimplotx, 'PrecondBandWidth', 0);
                end
                act2 = 1:length(paraAct);
                actN2 = act2;
                act2(paraAct==0) = [];  % Contains indices of variables to be optimized
                actN2(paraAct==1) = []; % Contains indices of variables not to be optimized

                truninitA = parainit(act2);
                truninitN = parainit(actN2); % Model still needs both types of variables to be calculated
                trunLB = paraLB(act2);
                trunUB = paraUB(act2);

                len = length(parainit);
                Id = eye(len);
                IdA = Id(:, act2);  % Active columns
                IdN = Id(:, actN2); % Inactive columns

                % Error function is model.RES, needs lenX
                lenX = zeros(length(X),1);
                for j=1:length(X)
                    lenX(j) = length(X{j});
                end

                trunFun = @(params) obj.RES(IdA*params(:)+IdN*truninitN(:), X, D, ma, lenX, BC); % 'params' = modified 'truninitA'

                % Optimization
                trunresults = lsqnonlin(trunFun, truninitA, trunLB, trunUB, optopt);

                % Replace old parameters and rescale them
                paraOpt = IdA*trunresults + IdN*truninitN;
                paraOpt = obj.scaleM*paraOpt;
                paraLB = obj.scaleM*paraLB;
                paraUB = obj.scaleM*paraUB;
                parainit = paraOpt;
            end
            obj.opt = obj.vec2para(parainit);  % Reverse index magic
            obj.scaleM = eye(length(parainit));
            obj.scaleMi = obj.scaleM;
            comptime = toc;

            if ~isempty(spath) % Save the local variables and data in workspace to file, if so desired.
                save([spath filesep 'res_after_calc.mat'])
            end
        end

        %% Converting parameter structs into vectors (index magic)
        % Vectors have the form (amt_1,1 amt_2,1, ... amt_n,m center_1,1 ..., lambda_n,m),
        % where 'xyz_i,j' is the value of the 'i'-th parameter 'xyz' in layer 'j'
        function vec = para2vec(obj, stru)
            lenC = numel(stru.center);
            vec = zeros(length(obj.ParaName)*lenC, 1);
            for k = 1:length(obj.ParaName)
                eval(['vec(1+(k-1)*lenC:k*lenC) = stru.' obj.ParaName{k} '(:);'])
            end
        end

        %% Converting vector into parameter structs (reverse index magic)
        function stru = vec2para(obj, vec) %#ok<STOUT>
            BigParaMatrix = reshape(vec, obj.lenPeak, length(obj.Ts), length(obj.ParaName)); %#ok<NASGU>
            for k = 1:length(obj.ParaName)
                eval(['stru.' obj.ParaName{k} '=BigParaMatrix(:,:,k);'])
            end
        end

        %% Auxiliary function for showing parameters
        function  pv = para2vecSHOW(obj)
            tName       = cell(obj.lenPeak*length(obj.ParaName)*length(obj.Ts),1);
            tParaLB     = zeros(obj.lenPeak*length(obj.ParaName)*length(obj.Ts),1);
            tParaINIT   = tParaLB;
            tParaUB     = tParaLB;
            tActive     = tParaLB;
            tParaOPT    = tParaLB;
            paraanz     = tParaLB;

            pa_idx = 1;
            l = obj.lenPeak;
            p = length(obj.ParaName);
            ind = 1;
            for i = 1:l
                for k = 1:p
                    eval(['cinit = obj.init.' obj.ParaName{k} ';'])
                    eval(['clower = obj.lower.' obj.ParaName{k} ';'])
                    eval(['cupper = obj.upper.' obj.ParaName{k} ';'])
                    eval(['cact = obj.act.' obj.ParaName{k} ';'])
                    eval(['copt = obj.opt.' obj.ParaName{k} ';'])

                    for j = 1:length(obj.Ts)

                        tName{ind}      = ['Mult' num2str(i) '_' obj.ParaName{k} '_' num2str(j)];
                        tParaLB(ind)    = clower(i,j);
                        tParaINIT(ind)  = cinit(i,j);
                        tParaUB(ind)    = cupper(i,j);
                        tActive(ind)    = logical(cact(i,j));
                        tParaOPT(ind)   = copt(i,j);
                        paraanz(ind)    = pa_idx;

                        ind = ind+1;
                    end
                    pa_idx = pa_idx+1;
                end
            end
            pv.Name     = tName(:);
            pv.ParaLB   = tParaLB(:);
            pv.ParaINIT = tParaINIT(:);
            pv.ParaUB   = tParaUB(:);
            pv.Active   = logical(tActive(:));
            pv.ParaOPT  = tParaOPT(:);
            pv.paraanz  = paraanz(:);
        end

        %% Showing and manipulating parameters
        function obj = showpara(obj)
            tf = figure('visible', 'off');

            % Create figure
            uf = uifigure;
            uf.Name=['Parameter Overview (Close to proceed)']; %#ok<NBRAK>
            uf.Position = [100 100 700 500];

            % Create table
            t=table('Size',[obj.lenPara 6],...
                'VariableTypes',{'string','double','double','double','logical','double'},...
                'VariableNames',{'Name','ParaLB','ParaINIT','ParaUB','Active','ParaOPT'});

            % Transform para matrices to single vector
            pv=obj.para2vecSHOW;
            t.Name      =pv.Name;
            t.ParaLB    =pv.ParaLB;
            t.ParaINIT  =pv.ParaINIT;
            t.ParaUB    =pv.ParaUB;
            t.Active    =pv.Active;
            t.ParaOPT   =pv.ParaOPT;
            paraanz     =pv.paraanz;

            % Put table in figure
            uit = uitable(uf);
            uit.Position=[10 10 680 480];
            uit.Data = t;
            uit.ColumnEditable = [false true true true true true];
            uf.CloseRequestFcn=@(src,call) my_closereq(src,call,uf,uit,tf); % Function is in \bin folder

            % Set cell edit function
            % Special cell edit to enable blockwise (de)activation of active paramters
            set(uit,'CellEditCallback',{@celledit,uf,paraanz});

            % Wait for figure to be closed
            uiwait(uf);

            % Write table data to output
            t=getappdata(tf,'tablecontent'); %#ok<NASGU>
            delete(tf);

            % Write changed parameters into model
            idx=1;
            for i=1:obj.lenPeak
                for kk=1:length(obj.ParaName)
                    for j=1:length(obj.Ts)
                        eval(['obj.init.' obj.ParaName{kk} '(' num2str(i) ',' num2str(j) ')=t.ParaINIT(' num2str(idx) ');'])
                        eval(['obj.lower.' obj.ParaName{kk} '(' num2str(i) ',' num2str(j) ')=t.ParaLB(' num2str(idx) ');'])
                        eval(['obj.upper.' obj.ParaName{kk} '(' num2str(i) ',' num2str(j) ')=t.ParaUB(' num2str(idx) ');'])
                        eval(['obj.act.' obj.ParaName{kk} '(' num2str(i) ',' num2str(j) ')=  t.Active(' num2str(idx) ');'])
                        eval(['obj.opt.' obj.ParaName{kk} '(' num2str(i) ',' num2str(j) ')=  t.ParaOPT(' num2str(idx) ');'])
                        idx=idx+1;
                    end
                end
            end
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variants (all) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Optimizing the 'center' and 'dist' values inside the reduced model
        function comptime = optimizeCenters(obj, opt, X, CP)
            tic;
            centers = obj.init.center(:); % Vector of length lenPeak * length(Ts)
            amts = obj.init.amt(:,1); % Same everywhere
            dists = obj.init.dist(:); % We need to optimize active dists as well
            actind = obj.act.dist(:);
            % Everything needs to be in vector-form for lsqnonlin
            LB_1 = min(X{1})*ones(size(centers));
            UB_1 = max(X{1})*ones(size(centers));
            LB_2 = 0.01*ones(size(dists(actind)));
            UB_2 = opt.L1opt.maxdist*ones(size(dists(actind)));
            LB = [LB_1; LB_2];
            UB = [UB_1; UB_2];


            nr_runs = length(opt.scheme);
            for i=1:nr_runs
                disp(['RUN ' num2str(i)]);

                switch opt.precon
                    case 1
                        optopt = optimoptions(@lsqnonlin, 'Display', 'Iter', 'MaxIter', opt.Init.MaxIter, 'MaxFunEvals', opt.Init.MaxFunEvals, 'TolFun', opt.Init.TolFun, 'TolX', opt.Init.TolX, 'plotfcn', @optimplotx);
                    case 2
                        optopt = optimoptions(@lsqnonlin, 'Display', 'Iter', 'MaxIter', opt.Init.MaxIter, 'MaxFunEvals', opt.Init.MaxFunEvals, 'TolFun', opt.Init.TolFun, 'TolX', opt.Init.TolX, 'plotfcn', @optimplotx, 'PrecondBandWidth', 0);
                end

                % If error function is model.CenterRES, it needs lenX
                lenX = zeros(length(X),1);
                for j=1:length(X)
                    lenX(j) = length(X{j});
                end

                centFun3 = @(params) obj.CenterRES3(params, amts, X, CP); % Third time's the charm.

                % Optimization, replaces old parameters
                params = lsqnonlin(centFun3, [centers; dists(actind)], LB, UB, optopt);

            end
            longcenters = params(1:obj.lenPeak*length(obj.Ts));
            longdists = params(obj.lenPeak*length(obj.Ts)+1:end);
            longerdists = zeros(size(longcenters));
            obj.opt.center = reshape(longcenters, obj.lenPeak, length(obj.Ts));
            ind = 1;
            for i = 1:length(actind)
                if actind(i)
                    longerdists(i) = longdists(ind);
                    ind = ind +1;
                else
                    longerdists(i) = -1;
                end
            end
            obj.opt.dist = reshape(longerdists, obj.lenPeak, length(obj.Ts));
            comptime = toc;
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variant 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Getting reduced and bounded model of this existing model (few timelayers & multiplets)
        function boundmod = getBoundedModel(obj, redTsidxIdx, multidx)
            redTsidx = obj.Tsidx(redTsidxIdx);

            boundmod = GLmodel(obj.Tfull(min(redTsidx):max(redTsidx)));

            boundmod.Ts = obj.Tfull(redTsidx);
            boundmod.Tsidx = redTsidx - min(redTsidx)+1; % In the reduced model the timelayers start at 'Tfull(min(redTsidx))'
            lenTs = length(boundmod.Ts);
            boundmod.lenPeak = length(multidx);
            boundmod.lenPara = length(boundmod.ParaName)*boundmod.lenPeak*lenTs;

            for k=1:length(boundmod.ParaName)   % Eval prevents writing 6 blocks of very similar code.
                eval(['boundmod.init.' obj.ParaName{k} '=obj.init.' obj.ParaName{k} '(multidx, redTsidxIdx);'])
                eval(['boundmod.lower.' obj.ParaName{k} '=obj.lower.' obj.ParaName{k} '(multidx, redTsidxIdx);'])
                eval(['boundmod.upper.' obj.ParaName{k} '=obj.upper.' obj.ParaName{k} '(multidx, redTsidxIdx);'])
                eval(['boundmod.act.' obj.ParaName{k} '=obj.act.' obj.ParaName{k} '(multidx, redTsidxIdx);'])
                eval(['boundmod.opt.' obj.ParaName{k} '=obj.opt.' obj.ParaName{k} '(multidx, redTsidxIdx);'])
            end

            for i=1:length(obj.constraints)
                boundmod.addConstraint(obj.constraints{i}.type, obj.constraints{i}.weight, obj.constraints{i}.var, obj.constraints{i}.Pidx);
            end
        end

        %% Inserting parameters from bounded model into big model
        function insertBoundedPara(obj,inp_opt,redTsidxIdx, multidx)  %#ok<INUSD>
            for k=1:length(obj.ParaName)
                eval(['obj.opt.' obj.ParaName{k} '(multidx ,redTsidxIdx)=inp_opt.' obj.ParaName{k} ';'])
            end
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variant 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Optimizing the model's first layer
        function obj = analyzeSingLayer(obj, Xinit, Dinit, specidx, options, BC)
            % Centers, dists and amts are already found.
            tidx = find(obj.Tsidx == specidx);
            centers = obj.opt.center(:,tidx);
            dists = obj.opt.dist(:, tidx);
            amts = obj.opt.amt(:,tidx);


            [scent, ord] = sort(centers);
            sdist = dists(ord);
            samts = amts(ord);

            L1opt = options.L1opt;
            bw = Xinit(1+options.bwdth)-Xinit(1);
            adapt = Xinit(1+options.L1opt.adapt)-Xinit(1);

            lefts = (Xinit(1)-1)*ones(length(scent),1); % Stores absolute values of local optimization boundaries, no indices here
            rights = zeros(length(scent),1);

            for j = 1:obj.lenPeak % Multiplets are easy
                if samts(j) ~= 1
                    lefts(j) = max(Xinit(1),scent(j) - max(bw, sdist(j))); % To the left of the whole multiplet
                    rights(j) = min(Xinit(end), scent(j) + max(bw*samts(j), sdist(j)*samts(j))); % To the right
                end
            end

            if (obj.lenPeak == 1 && samts(1) == 1) % A singular peak is easy
                lefts(1) = max(Xinit(1), scent(j) - adapt);
                rights(1) = min(Xinit(end), scent(j) + adapt);
            else
                for j=find(lefts < Xinit(1))' % Remaining singlets:
                    if j == 1
                        lefts(j) = max(Xinit(1), scent(j)-adapt); % Don't overshoot, and on right side dont overshoot 'adapt' or next peak as well
                        rights(j)= min([Xinit(end), scent(j)+adapt, (scent(j)+scent(j+1))/2]);
                    elseif j == obj.lenPeak
                        if rights(j-1) < scent(j) % Don't overshoot 'rights(j-1)' as well, since it denotes the end of a a multiplet
                            lefts(j) = max([Xinit(1), scent(j)-adapt, rights(j-1)]);
                            rights(j) = min(Xinit(end), scent(j)+adapt);
                        else  % Is only the case if the current peak is inside of another multiplet
                            lefts(j) = max(Xinit(1), scent(j)- bw); % We should have 'bw' << 'adapt'
                            rights(j) = min(Xinit(end), scent(j)+bw);
                        end

                    else
                        if rights(j-1) < scent(j)
                            lefts(j) = max([Xinit(1), scent(j)-adapt, rights(j-1)]);
                            rights(j)= min([Xinit(end), scent(j)+adapt, (scent(j)+scent(j+1))/2]);
                        else
                            lefts(j) = max(Xinit(1), scent(j)- bw);
                            rights(j) = min(Xinit(end), scent(j)+bw);
                        end

                    end
                end
            end


            optopt = optimset('Display','iter','MaxIter',options.Init.MaxIter,'MaxFunEvals',options.Init.MaxFunEvals,'TolFun',options.Init.TolFun,'TolX',options.Init.TolX);

            % Optimize 'dist', 'hw', 'height' and 'lambda' at once, first for each peak individually
            distinit = zeros(obj.lenPeak, 1);
            C = zeros(obj.lenPeak,1);
            for i=1:obj.lenPeak
                distinit(i) = sdist(i); % Set for all peaks, only used for true multiplets
                cidx = find(Xinit >= scent(i), 1);
                C(i) = 0.8*Dinit(cidx)*nchoosek(samts(i)-1, ceil((samts(i)-1)/2));
            end
            ainit = 0.01;
            lambinit = 0.1;
            LB_1 = [0,0,L1opt.minlamb];  % If amt == 1
            UB_1 = [L1opt.maxhw, inf, L1opt.maxlamb];
            LB_2 = [0,0,0,L1opt.minlamb]; % Else
            UB_2 = [L1opt.maxdist, L1opt.maxhw, inf, L1opt.maxlamb];

            for i = 1:obj.lenPeak
                lb = find(Xinit >= lefts(i), 1);
                rb = find(Xinit >= rights(i), 1);
                amt = samts(i);
                center = scent(i);
                if amt == 1
                    smol_fun = @ (para) (Dinit(lb:rb) -glmult(Xinit(lb:rb), amt, center, -1, para(1), para(2), para(3), BC)); % Errorfunction is difference from data to model
                    [para_opt, ~] = lsqnonlin(smol_fun, [ainit, C(i), lambinit], LB_1, UB_1, optopt);

                    % Write values into correct slot (respecting the order)
                    obj.opt.hw(ord(i),tidx) = para_opt(1);
                    obj.opt.height(ord(i), tidx) = para_opt(2);
                    obj.opt.lambda(ord(i), tidx) = para_opt(3);
                else
                    smol_fun = @ (para) (Dinit(lb:rb) -glmult(Xinit(lb:rb), amt, center, para(1), para(2), para(3), para(4), BC)); % Errorfunction is difference from data to model
                    [para_opt, ~] = lsqnonlin(smol_fun, [distinit(i), ainit, C(i), lambinit], LB_2, UB_2, optopt);

                    obj.opt.dist(ord(i), tidx) = para_opt(1); % Otherwise there's the optimized value
                    obj.opt.hw(ord(i), tidx) = para_opt(2);
                    obj.opt.height(ord(i), tidx) = para_opt(3);
                    obj.opt.lambda(ord(i), tidx) = para_opt(4);
                end
            end

            clower.amt       = obj.opt.amt(:,tidx);
            cupper.amt       = obj.opt.amt(:,tidx);
            cact.amt         = false(size(cupper.amt)); % Don't optimize 'amt'.

            clower.center    = zeros(size(obj.opt.center(:,tidx))); % Only temporary, no worries
            cupper.center    = zeros(size(obj.opt.center(:,tidx)));
            cact.center      = true(size(cupper.center));

            clower.dist      = zeros(size(obj.opt.dist(:,tidx)));
            cupper.dist      = options.L1opt.maxdist*ones(size(obj.opt.dist(:,tidx)));
            cact.dist        = true(size(cupper.dist)); % Only temporarily fully active

            clower.hw        = zeros(size(obj.opt.hw(:,tidx)));
            cupper.hw        = L1opt.maxhw*ones(size(obj.opt.hw(:,tidx)));
            cact.hw          = true(size(cupper.hw));

            clower.height    = zeros(size(obj.opt.height(:,tidx)));
            cupper.height    = inf*ones(size(obj.opt.height(:,tidx)));
            cact.height      = true(size(cupper.height));

            clower.lambda    = L1opt.minlamb*ones(size(obj.opt.lambda(:,tidx)));
            cupper.lambda    = L1opt.maxlamb*ones(size(obj.opt.lambda(:,tidx)));
            cact.lambda      = true(size(obj.opt.lambda(:,tidx)));

            for i=1:obj.lenPeak % See, i told you not to worry
                clower.center(ord(i)) = max(scent(i)-Xinit(1+L1opt.maxcenter)-Xinit(1), Xinit(1));
                cupper.center(ord(i)) = min(scent(i)+Xinit(1+L1opt.maxcenter)-Xinit(1), Xinit(end));
                if samts(i) == 1
                    cact.dist(ord(i)) = false; % Again, this works
                end
            end
            % Bring parameters in single-layer-form for model
            for param = obj.ParaName
                eval(['inits.' param{1} '=obj.opt.' param{1} '(:, tidx);'])
            end
            model = GLmodel(obj.Tfull(1), inits, clower, cupper, cact, 1); % Creates auxiliary model of a single time layer to optimize
            model.optimizeModel(options, {Xinit}, {Dinit}, max(Dinit), [], BC);
            % Bring parameters back out of single-layer-form
            for param = obj.ParaName
                eval(['obj.opt.' param{1} '(:,tidx) = model.opt.' param{1} ';'])
            end

            % ~ Nice plots x7
            if options.info >= 3
                figure;
                plot(Xinit, Dinit, 'k');
                set(gca, 'XDir', 'reverse');
                hold on
                for i=1:obj.lenPeak
                    plot(Xinit, glmult(Xinit, initLayer.amt(i), initLayer.center(i), initLayer.dist(i), initLayer.hw(i), initLayer.height(i), initLayer.lambda(i), BC), 'r--')
                end
            end
        end
    end
end