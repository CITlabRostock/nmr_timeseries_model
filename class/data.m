classdef data < handle
    % Stores data of the time series
    properties
        location    % Folder path - location of dat-files
        name        % Data name - no computational meaning
        comment     % Data comments - no computational meaning
        clim        % Scaling factor of colour gradient - no computational meaning
        D           % Cell of vectors - D{j}(i) = NMR data for T(j) and X{j}(i)
        T           % Vector of timelayers
        X           % Cell of vectors - X{j}(i) = frequncy for T(j) and X{j}(i)
        ma          % Maximum of every spectrum
        scal        % Scaling factor Draw = scal*D, 'D' is scaled to 1
        CP          % Cell of vectors similar to D{j}(i) containing probabilities for center optimization
        timenoisy   % Boolean that determines whether data needs smoothing through time dimension
        
    end
    
    methods
        %% Constructor 
        % Assumptions: 'D' contains the Data; either the rows represent each spectrum or the columns do
        % 'x' or 'X' contains the measured points in frequency domain
        % 't' or 'T' contains the measured points in time
        function obj = data(file,name,clim,comment, tn)
            
            obj.location=file;
            if exist('clim', 'var')
                obj.clim = clim;
            else
                obj.clim = [0,1];
            end
            if exist('comment', 'var')
                obj.comment = comment;
            else
                obj.comment = '';
            end
            if exist('name','var')
                obj.name=name;
            else
                obj.name = 'noname';
            end
            if exist('tn','var')
                obj.timenoisy = tn;
            else
                obj.timenoisy = 0;
            end
            % Loads specified variables into workspace
            load(file,'D','t','x','T','X');
            try
                t=T;
            catch
            end
            try
                x=X;
            catch
            end
                                        
            if size(D,1)>size(D,2) % We want to have #rows=#spectra, #columns=#freq. points. Usually, #spectra << #freq. points
                D = D';
            end
            
            if length(t) > size(D,1) % Might happen if you don't SINC all time layers
                t = t(1:size(D,1));
            end

            % Start the data transcription
            obj.scal = max(max(D)); 
            obj.T = t;
            obj.X={};
            obj.D={};
            obj.ma = zeros(1,length(obj.T));
            [x, order] = sort(x); % 'x' might be unsorted or in wrong direction .
            for i=1:length(obj.T)
                obj.X{end+1} = x(:);    % One could leave 'X' one-dimensional, but alas.
                temp = D(i,order)/obj.scal;
                obj.D{end+1}=temp(:); % Also sort 'D' that way, if necessary.
                obj.ma(i)=max(abs(temp));
            end
            obj.CP = cell(size(obj.X));
        end
        
        %% Copying instances of 'data' without overwriting
        function copied_data = copy(obj)
            loc = obj.location;
            nam = obj.name;
            cli = obj.clim;
            cmnt = [obj.name, 'w/ model subtracted'];
            copied_data = data(loc,nam,cli,cmnt);
            copied_data.T = obj.T;
            copied_data.D = obj.D;
            copied_data.X = obj.X;
            copied_data.ma = obj.ma;
            copied_data.scal = obj.scal;
            copied_data.CP = obj.CP;
            copied_data.timenoisy = obj.timenoisy;
        end
        %% Cutting off edges of 'data' (time-wise, spectra Tmin<=t<=Tmax remain)
        function obj = restrictT(obj,Tmin,Tmax)
            idx = 1:length(obj.T);
            idx(max(obj.T<Tmin, obj.T>Tmax)) = [];
            obj.T = obj.T(idx);
            obj.D = obj.D(idx);
            obj.X = obj.X(idx);
            obj.ma = obj.ma(idx);
            obj.CP = obj.CP(idx);
        end
        
        %% Thinning out 'data' (time-wise, to percent%)
        function obj = thinoutT(obj,percent)
            idx=round(linspace(1,length(obj.T), ceil(percent/100*length(obj.T))));
            if(length(idx)< 5)
                warning(['thinoutT: new data set has too few spectra (<5). Aborting.']) %#ok<NBRAK>
                return
            end
            obj.T = obj.T(idx);
            obj.D = obj.D(idx);
            obj.X = obj.X(idx);
            obj.ma = obj.ma(idx);
            obj.CP = obj.CP(idx);
        end
        
        %% Cutting off edges of 'data' (frequency-wise, columns with Xmin<=x<=Xmax remain) 
        function obj = restrictX(obj,Xmin,Xmax)
            for i=1:length(obj.X)
                idx = 1:length(obj.X{i});
                idx(max(obj.X{i}<Xmin, obj.X{i}>Xmax)) = [];
                
                obj.D{i} = obj.D{i}(idx);
                obj.X{i} = obj.X{i}(idx);
                obj.ma(i) = max(abs(obj.D{i}));
                if ~isempty(obj.CP{i})
                    obj.CP{i} = obj.CP{i}(idx);
                end
            end
        end
        
        %% Thinning out 'data' (frequency-wise, to percent%)
        function obj = thinoutX(obj, percent)
           for i=1:length(obj.X)
               idx = round(linspace(1,length(obj.X{i}), ceil(percent/100*length(obj.X{i}))));
               obj.D{i} = obj.D{i}(idx);
               obj.X{i} = obj.X{i}(idx);
               obj.ma(i) = max(abs(obj.D{i}));
               if ~isempty(obj.CP{i})
                    obj.CP{i} = obj.CP{i}(idx);
               end
           end
        end
        
        %% Printing available information 
        function print(obj,indent) % Larger 'indent' results in more space before information is printed
            if nargin==1
                indent=0;
            end
            data=obj;
            empStr(1:indent)=' ';
            disp([empStr '              Data: ' data.name])
            disp([empStr '          Comments: ' data.comment])
            disp([empStr '          Location: ' data.location])
            disp([empStr '  Length time grid: ' num2str(length(data.T))])
            try
                disp([empStr 'Length CShift grid: ' num2str(length(data.X{1}))])
            catch
            end
        end
        
        %% Plotting 'data'
        function plot(obj, xlim, ylim)
            col1=[0    0.4470    0.7410];       % colour scheme
            col2=[0.8500    0.3250    0.0980];
            
            figure;
            try
                xmin = find(obj.X{1}>=xlim(1), 1,  'first');
                xmax = find(obj.X{1}<=xlim(2), 1, 'last');
                
            catch
                xmin = 1;
                xmax = length(obj.X{1});
            end
            try
                ymin = find(obj.T>=ylim(1),1,'first');
                ymax = find(obj.T<=ylim(2),1,'last');
            catch
                ymin = 1;
                ymax = length(obj.T);
            end
            try
                Dm=zeros(xmax-xmin+1,ymax-ymin+1);
                for j=ymin:ymax
                    Dm(:,j-ymin+1)=obj.D{j}(xmin:xmax);
                end
                
                surf(obj.T(ymin:ymax),obj.X{1}(xmin:xmax),Dm)
                shading interp
                view([-60 45]);
            catch
                warning('Data not rectangular')
                for j=ymin:ymax
                    l=(j-ymin)/(ymax-ymin);
                    plot3(obj.X{j}(xmin:xmax),obj.T(j)*ones(xmax-xmin+1,1),obj.D{j}(xmin:xmax),'color',l*col2+(1-l)*col1);
                    hold on
                end
            end
            axis tight
        end
        
        %% Checking data
        function reg = check(obj)
            b=1;
            irreg = 0;
            irrec = 0;
            for i=1:length(obj.T)-1
                if length(obj.X{i})==length(obj.X{i+1})
                    if max(abs(obj.X{i}-obj.X{i+1}))~=0
                        b=0;
                        irreg = i;
                    end
                else
                    b=0;
                    irrec = i;
                end
            end
            if b
                disp('Data check: Data is regular and rectangular')
                reg = 1;
            else
                disp('Data check: Data is not regular')
                if irreg == 0
                    disp(['Lengths of X{' num2str(irrec) '} and X{' num2str(irrec+1) '} differ']) 
                else
                    disp(['Point values of X{'  num2str(irreg) '} and X{' num2str(irreg+1) '} differ'])
                end
                reg = 0;
            end
        end
        
        %% Calculating the sum of all data points, relevant for model quality metrics
        function mass = calcsum(obj)
            mass = 0;
            for j = 1:length(obj.D)
                l = length(obj.D{j});
                for i=1:l
                    mass = mass+max(0,obj.D{j}(i));
                end
            end
        end
        %% Smoothing data over time
        function obj = smoothTime(obj, tbw)
            dm = zeros(length(obj.D),length(obj.D{1}));
            for i=1:length(obj.D)
                dm(i,:)= obj.D{i};
            end
            dm_smooth = smoothdata(dm,1,'movmean',tbw); % Works columnwise, pretty neat
            
            for j=1:length(obj.D)
                obj.D{j} = dm_smooth(j,:)'; % Careful with rows and columns
            end
        end
        
        %% Smoothing data in the rough direction of peak movements
        function obj = smoothTimeDetailed(obj, tbw, offpx, bds)
            % Reading the data
            dm = zeros(length(obj.D),length(obj.D{1}));
            dm_smooth = dm;
            for i=1:length(obj.D)
                dm(i,:)= obj.D{i};
            end
            
            % Smoothing procedure 
            numtl = size(dm,1);
            numpx = size(dm,2);
            for i=1:numtl
                lower = min(floor(tbw/2), i-1);
                upper = min(floor(tbw/2), numtl-i);
                for j=1:numpx
                    region = max(1,find(bds(i,:) >= j,1)-1); % Special case j = 1 is covered by the maximum
                    realoffset = offpx(i-lower:i+upper, region)+j-offpx(i, region);
                    realoffset = max(realoffset, bds(i-lower:i+upper, region)); % Cutting off too low values at bounds
                    realoffset = min(realoffset, bds(i-lower:i+upper, region+1)); % Cutting off too high values at bounds
                    
                    smoothsum = 0;
                    for k = 1:length(realoffset)
                        smoothsum = smoothsum + dm(i-lower+k-1, realoffset(k)); % Moving average Filter, can be adapted.
                    end
                    dm_smooth(i,j) = smoothsum/(lower+upper+1);
                end
            end
            
            % Change data
            for j=1:length(obj.D)
                obj.D{j} = dm_smooth(j,:)'; % Careful with rows and columns
            end
        end
        
        %% Smoothing data near the peaks w/ consideration of peak movement to rid the errors due to crossing solvent peaks
        function obj = smoothTimeSausage(obj, tbw, swdth, cpx) 
            % Reading the data
            dm = zeros(length(obj.D),length(obj.D{1}));
            numruns = dm; % # Smoothing runs per pixel, mostly zero in which case dm is copied
            dm_smooth = dm; % Allocate correct dimensional matrix first
            for i=1:length(obj.D)
                dm(i,:)= obj.D{i};
            end
            
            % Smoothing procedure
            numtl = size(cpx,1);
            numpk = size(cpx,2);
            numpx = size(dm,2);
            shw = floor(swdth/2);
            
            % Preparing the Smoothing to accomodate multiple smoothing runs on single pixels
            for i=1:numtl
                for j=1:numpk % Need to do that for each peak
                    left = max(1, cpx(i,j)-shw);
                    right = min(numpx, cpx(i,j)+shw); 
                    for k = left:right % Wobbly rectangle aka sausage
                        numruns(i,k) = numruns(i,k) + 1;
                    end
                end
            end
            dm_smooth(numruns == 0) = dm(numruns == 0); % Smooth data where nothing is smoothed is just the data
            
            % Begin actual smoothing calculations
            for i=1:numtl
                lower = min(floor(tbw/2), i-1);
                upper = min(floor(tbw/2), numtl-i);
                for j=1:numpk % Need to do that for each peak
                    left = max(1, cpx(i,j)-shw);
                    right = min(numpx, cpx(i,j)+shw); 
                    for k = left:right % Wobbly rectangle, aka sausage
                        smoothsum = 0;
                        for l = (i-lower):(i+upper) % Window in time
                            pxpos = k + cpx(l,j)-cpx(i,j); % Peak offset in the time-window
                            pxpos = max(1,min(numpx,pxpos));
                            smoothsum = smoothsum + dm(l, pxpos);
                        end
                        if numruns(i,k) > 1
                            ratio = 0;
                            for m = 1:numruns(i,k) % Change m to proper number of the peak and numruns to 3D-matrix
                                ratio = ratio + dm(i-lower, max(1,min(numpx,k+cpx(i-lower,m)-cpx(i,m))))/2 + dm(i+upper, max(1,min(numpx,k+cpx(i+upper,m)-cpx(i,m))))/2;
                            end
                        ratio = (dm(i-lower, max(1,min(numpx, k+cpx(i-lower,j)-cpx(i,j))))/2+dm(i+upper, max(1,min(numpx,k+cpx(i+upper,j)-cpx(i,j))))/2)/ratio;
                        else
                            ratio = 1;
                        end
                        dm_smooth(i,k) = dm_smooth(i,k) + (smoothsum/(lower+upper+1))*ratio;
                    end
                end
            end
            
            % Change data
            for j=1:length(obj.D)
                obj.D{j} = dm_smooth(j,:)'; % Careful with rows and columns
            end
        end
        
        %% Smoothing data near the peaks w/ consideration of peak movement to rid the errors due to crossing solvent peaks
        function obj = smoothTimeSausage2(obj, tbw, swdth, d, cpx)
            if mod(tbw, 2) == 0
                warning('tbwdth must be odd. Abort smoothing')
                return
            end
            % Reading in from data
            dm = zeros(length(obj.D),length(obj.D{1}));
            for i=1:length(obj.D)
                dm(i,:)= obj.D{i};
            end
            dm_smooth = dm;
                        
            % Smoothing procedure
            numtl = size(cpx,1);
            numpk = size(cpx,2);
            numpx = size(dm,2);
            shw = floor(swdth/2);
            
            % Begin actual smoothing calculations
            for i=1:numpk
                for shiftind = -shw:0 % Splitting of the loops for overwriting superfluous edgeinformation
                    inds = max(1, cpx(:,i)+shiftind);
                    skewcol = obj.calcskewcol(inds);
                    skewcol2 = mysgfilt(d,tbw, skewcol);
                    for t=1:numtl
                        dm_smooth(t,inds(t)) = skewcol2(t);
                    end
                end
                for shiftind = shw:-1:1
                    inds = min(numpx, cpx(:,i)+shiftind);
                    skewcol = obj.calcskewcol(inds);
                    skewcol2 = mysgfilt(d,tbw, skewcol);
                    masmooth = smoothdata(skewcol, 'movmean', tbw);
                    thw = (tbw-1)/2;
                    skewcol2([1:thw, end-thw+1:end]) = masmooth([1:thw, end-thw+1:end]); % Edges via moving average, innards via Savitzky-Golay
                    for t=1:numtl
                        dm_smooth(t,inds(t)) = skewcol2(t);
                    end
                end
            end
            
            % Change data
            for j=1:length(obj.D)
                obj.D{j} = dm_smooth(j,:)'; % Careful with rows and columns
            end
        end
        
        %% Smoothing data via multiplying factor to whole spectrum, UNUSED
        function obj = smoothTimeFactor(obj, facs)
            for i=1:length(obj.T)
                obj.D{i} = obj.D{i}*facs(i);
            end
        end
        
        %% Calculating the skewed columns
        function vec = calcskewcol(obj,inds)
            if length(inds) == 1 && length(obj.T) > 1
                inds = inds*ones(length(obj.T),1);
            end
            vec = zeros(length(obj.T),1);
            for t=1:length(obj.T)
                vec(t) = obj.D{t}(inds(t));
            end
        end
    end
end
