classdef (Abstract) model < handle
    % Superclass for all models. Instances of Subclasses contain all
    % parameters of modeled timelayers. Model calculation is done by
    % subclasses. Superclass calculates residuals of
    % (model-data)^2 + violation of constraint errors.

    properties (Abstract)
        mname                       % Predefined name of model - don't change
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

        ParaName                    % Names of parameters - don't change

    end

    properties
        constraints={};             % Cell of constraints
    end

    methods (Abstract)
        MOD(obj,para,RDobj)         % Obtains modeled values in 'D{i}(j)'-format
        plot(obj,RDobj)             % Plots model (no surprise)
    end

    methods
        %% Calculate residuals to input data
        function R = RES(obj, para, X, D, ma, lenX, BC)
            M = obj.MOD(para, X, BC);     % 'M' has same format as 'D', 'X'

            R = zeros(sum(lenX),1);   % lenX(i) = length(X{i})

            % Residuals from the model
            ind = 0;
            isayso = 0; % Testing purposes only
            for i = 1:length(D)
                if isayso
                    vec = 1/ma(i)*(D{i}-M{i}); %#ok<*UNRCH>
                    for j=1:length(vec)
                        if vec(j) >= 0
                            vec(j) = sqrt(vec(j));
                        else
                            vec(j) = -sqrt(-vec(j));
                        end
                    end
                else
                    vec = 1/ma(i)*(D{i}-M{i});
                end
                R(ind+1:ind+lenX(i)) = vec; % no absoulte residuals because lsqnonlin demands it
                ind = ind + lenX(i);

            end
            % Residuals from the constraints
            for i=1:length(obj.constraints)
                c=obj.constraints{i}.eval(obj.scaleM*para);
                R = [R; c(:)]; %#ok<AGROW>
            end
        end

        %% Calculate residuals for center optimization, UNUSED VARIANT
        function R = CenterRES(obj, centers, amts, dists, X, CP, lenX)
            R = zeros(sum(lenX),1);
            centerfull = zeros(obj.lenPeak, length(obj.Tfull));
            if length(obj.Ts) == 1
                centerfull(:,1) = centers;
            else
                for i=1:obj.lenPeak
                    centerfull(i,:) =  pchip(obj.Ts, centers(i:obj.lenPeak:end), obj.Tfull);
                end
            end
            M = cell(size(CP));
            for i=1:length(CP)
                peaks = zeros(1,sum(amts)); % Breaking up multiplets
                pind = 1;
                for k=1:length(amts)
                    peaks(pind) = centerfull(k,i);
                    pind = pind+1;
                    for j = 2:amts(k)
                        peaks(pind) = centerfull(k,i)+(j-1)*dists(k);
                        pind = pind + 1;
                    end
                end
                peakidx = zeros(size(peaks));
                for k = 1:length(peaks)
                    peakidx(k) = find(X{i} >= peaks(k),1);
                end
                temp = zeros(size(X{i}));
                temp(peakidx) = 1;
                M{i} = smoothe_fct(temp); % Smoothes the discrete metric for better gradient optimization
            end
            ind = 0;
            for i = 1:length(X)
                R(ind+1:ind+lenX(i)) = CP{i}-M{i}; % Again, no absolute residuals
                ind = ind + lenX(i);
            end
        end

        %% Calculating residuals for center optimization, UNUSED VARIANT
        function R = CenterRES2(obj, centers, amts, dists, X, CP)
            R = zeros(sum(amts)*length(X),1); % #peaks x #timelayers
            centerfull = zeros(obj.lenPeak, length(obj.Tfull));
            if length(obj.Ts) == 1
                centerfull(:,1) = centers;
            else
                for i=1:obj.lenPeak
                    centerfull(i,:) =  pchip(obj.Ts, centers(i:obj.lenPeak:end), obj.Tfull);
                end
            end
            for i=1:length(CP)
                peaks = zeros(1,sum(amts)); % Breaking up multiplets
                pind = 1;
                for k=1:length(amts)
                    peaks(pind) = centerfull(k,i);
                    pind = pind+1;
                    for j = 2:amts(k)
                        peaks(pind) = centerfull(k,i)+(j-1)*dists(k);
                        pind = pind + 1;
                    end
                end
                for j=1:length(peaks)
                    R((i-1)*length(peaks)+j) = min(abs(CP{i}-peaks(j))); % Index magic. Could be more efficient.
                end
            end
            for i=1:length(obj.constraints)
                c=obj.constraints{i}.eval(obj.scaleM*para);
                R = [R; c(:)]; %#ok<AGROW>
            end
        end

        %% Calculating residuals for Center optimization, IN USE
        function R = CenterRES3(obj, params, amts, X, CP)
            R = zeros(sum(amts)*length(X),1); % #peaks x #timelayers
            centerlen = length(obj.Ts)*obj.lenPeak;
            distlen = length(params)-centerlen;
            centerfull = zeros(obj.lenPeak, length(obj.Tfull));
            distfull = -1*ones(obj.lenPeak, length(obj.Tfull));
            centers = params(1:centerlen);
            dists = params(centerlen+1:end);
            if length(obj.Ts) == 1
                centerfull(:,1) = centers;
                distfull(amts > 1,1) = dists;
            else
                ind = 1;
                for i=1:obj.lenPeak
                    centerfull(i,:) =  pchip(obj.Ts, centers(i:obj.lenPeak:end), obj.Tfull);
                    if amts(i) ~= 1
                        distfull(i,:) = pchip(obj.Ts, dists(ind:distlen/length(obj.Ts):end), obj.Tfull); % Retrieve information from dists
                    else
                        distfull(i,:) = -1*ones(size(distfull(i,:))); % Don't retrieve information from dists
                    end
                end
            end
            for i=1:length(CP)
                peaks = zeros(1,sum(amts)); % Breaking up multiplets
                pind = 1;
                for k=1:length(amts)
                    peaks(pind) = centerfull(k,i);
                    pind = pind+1;
                    for j = 2:amts(k)
                        peaks(pind) = centerfull(k,i)+(j-1)*distfull(k, i);
                        pind = pind + 1;
                    end
                end
                for j=1:length(peaks)
                    R((i-1)*length(peaks)+j) = min(abs(CP{i}-peaks(j))); % Index magic. Calculates the closest CP to any peak.
                end
            end
            for i=1:length(obj.constraints)
                c=obj.constraints{i}.eval(obj.scaleM*para);
                R = [R; c(:)]; %#ok<AGROW>
            end
        end

        %% Get information for creating constraints
        function st = getConstraintInfo(obj)
            st.name = obj.ParaName;
            st.range = {};
            lenSet = numel(obj.init.Center);
            for k = 1:length(obj.ParaName)
                st.range{end+1}= 1 + (k-1)*lenSet:kk*lenSet;
                % Allocates indices. range{j}(i) is the i-th peak's j-th parameter
            end
        end

        %%  Print information for existing constraints
        function printConstraintInfo(obj)
            constr = obj.constraints;
            disp(' ')
            disp(['Applied constraints for model ' obj.mname]);
            for i=1:length(constr)
                disp([num2str(i) ': ' constr{i}.var ' ' constr{i}.name ' ' constr{i}.weight ]);
            end
        end

        %% Add constraint for parameter (types: monotone/increase/decrease)
        function obj = addConstraint(obj,type,weight,var,Pidx)
            % disp('Add Constraint - START')
            % 'type' is one of monotone/increase/decrease
            % 'var' refers to amt/center/dist/hw/height/lambda
            if nargin ~= 5
                error('4 arguments needed - (1) constraint type, (2) weight, (3) {amt, center, dist, hw, height, lambda}, (4) PeakIndex ')
            end

            if weight<=0
                error('The weight has to be positive.')
            end

            b=0;
            for ch={'constant','increase','decrease'}
                if strcmp(ch,type)
                    b=1;
                    break
                end
            end

            if b==0
                error("Unknown constraint - options: 'constant', 'increase', 'decrease'")
            end


            st=obj.getConstraintInfo;
            % Check var
            b=0;
            for i=1:length(st.name)
                if strcmp(st.name{i},var)
                    stidx=i;
                    b=1;
                    break
                end
            end
            parastart=st.range{stidx}(1)-1;


            if b==0
                error('Third argument doesnt fit any of the model parameters.')
            end

            con.type=type;
            con.weight=weight;
            con.var=var;
            con.Pidx=Pidx; % Peakindex



            idxMat=reshape(1:obj.lenPeak*length(obj.Ts),obj.lenPeak,length(obj.Ts));


            switch type
                case 'constant'
                    con.name='constant';
                    con.eval=@(para) weight*min(...
                        abs(min(diff(para(parastart+idxMat(Pidx,:)))./para(parastart+idxMat(Pidx,1:end-1)),0)),...
                        abs(max(diff(para(parastart+idxMat(Pidx,:)))./para(parastart+idxMat(Pidx,1:end-1)),0))...
                        );

                case 'increase'
                    con.name='increase';
                    con.eval=@(para) weight*min(diff(para(parastart+idxMat(Pidx,:)))./para(parastart+idxMat(Pidx,1:end-1)),0);

                case 'decrease'
                    con.name='decrease';
                    con.eval=@(para) weight*max(diff(para(parastart+idxMat(Pidx,:)))./para(parastart+idxMat(Pidx,1:end-1)),0);


                otherwise
                    error("Unknown constraint - options: 'constant', 'increase', 'decrease'")
            end

            obj.constraints{end+1}=con;
            %disp('Add Constraint - END')

        end
    end
end

