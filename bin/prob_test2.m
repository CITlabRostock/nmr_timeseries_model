%% UNUSED predecessor of prob_test3.m 
function probs = prob_test2(d1, bw, adapt, nop, amp_min)
    probs = zeros(size(d1));
    ind = bw;
    if d1(bw) == 0
        d1(bw) = eps;
    end
    maxs = {};  % There always is a pair of mins and maxs. Additionally, the corresponting max always comes first
    mins = {};
    while ind <= length(probs) - bw
        currcase = sign(d1(ind));
        switch currcase % If a 0 occurs in 'd1', we shall stay in the current case until the sign strictly switches
            case 1
                maxs{end+1} = [d1(ind), ind]; %#ok<AGROW>
                while  ind<= length(probs) && sign(d1(ind)) >= 0 
                    if d1(ind) > maxs{end}(1)
                        maxs{end} = [d1(ind), ind];
                    end
                    ind = ind + 1;
                end
            case -1
                mins{end+1} = [d1(ind), ind]; %#ok<AGROW>
                while ind<= length(probs) && sign(d1(ind)) <= 0 
                    if d1(ind) < mins{end}(1)
                        mins{end} = [d1(ind), ind];
                    end
                    ind = ind +1;
                end
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
        scores(i) = maxs{i}(1)-mins{i}(1); % Adds both up ('mins' contains only negative values)
    end
    [scores, order] = sort(scores, 'desc'); % Largest scores first
    if length(scores) < nop % If there are not enough peaks, there will be fewer spikes in probability
        nop = length(scores);
    end
    peaks = zeros(nop, 1);
    for i = 1 : nop
        if scores(i) > amp_min
            peaks(i) = ceil((mins{order(i)}(2)+maxs{order(i)}(2))/2); % Middle between min and max
        end
    end
    peaks(peaks == 0) = [];
    probs(peaks) = 1;
    probs = smoothe_fct(probs, adapt);
end