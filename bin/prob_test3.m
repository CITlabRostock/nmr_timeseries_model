 % Instead of outputting a continuous function, we take the proposed peaks 
 % (in ppm, not indices) instead and get the minimal distance to the
 % proposed peaks
function probs = prob_test3(d0, d1, bw, nop, amp_min)
    ind = bw;
    if d1(bw) == 0
        d1(bw) = eps;
    end
    maxs = {};  % There always is a pair of mins and maxs. Additionally, the corresponting max always comes first
    mins = {};
    posnegind = {}; % Indices where we switch from >=0 to <0.
    while ind <= length(d1) - bw
        currcase = sign(d1(ind));
        switch currcase % If a 0 occurs in 'd1', we shall stay in the current case until the sign strictly switches
            case 1
                maxs{end+1} = [d1(ind), ind]; %#ok<*AGROW>
                while  ind<= length(d1) && sign(d1(ind)) >= 0 
                    if d1(ind) > maxs{end}(1)
                        maxs{end} = [d1(ind), ind];
                    end
                    ind = ind + 1;
                end
                posnegind{end+1} = ind-1; % If 'posnegind' grows when we finish >0, this entry will be ignored 
            case -1
                mins{end+1} = [d1(ind), ind];
                while ind<= length(d1) && sign(d1(ind)) <= 0 
                    if d1(ind) < mins{end}(1)
                        mins{end} = [d1(ind), ind];
                    end
                    ind = ind +1;
                end
            case 0  % Only occurs if first checked value of 'd1' is 0
                ind = ind +1;
        end
    end
    if mins{1}(2) < maxs{1}(2) % Repairing extrema, if necessary
        mins = mins(2:end);
    end
    if maxs{end}(2) > mins{end}(2)
        maxs = maxs(1:end-1);
    end
    scores = zeros(length(maxs),1);
    for i=1:length(scores)
        scores(i) = maxs{i}(1)-mins{i}(1); % Adds both up
    end
    [scores, order] = sort(scores, 'desc'); % Largest scores first
    if length(scores) < nop % If there are not enough peaks, we'll only find these peaks
        nop = length(scores);
    end
    peaks = zeros(nop, 1);
    for i = 1 : nop
        if scores(i) > amp_min
            peaks(i) = posnegind{order(i)};  % Switch index
        end
    end
    peaks(peaks == 0) = [];
    % Postprocessing of peaks
    for i=1:length(peaks)
        ind = peaks(i);
        while d0(ind)<d0(ind+1) && ind < peaks(i)+floor(bw/2)
            ind = ind+1;
        end
        while d0(ind)<d0(ind-1) && ind > peaks(i)-floor(bw/2)
            ind = ind-1;
        end
        peaks(i) = ind;
    end
    probs = peaks;
end