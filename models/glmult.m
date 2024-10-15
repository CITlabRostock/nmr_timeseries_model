function [res] = glmult(x, amt, center, dist, hw, height, lambda, BC)
% lambda == 0: Lorentzian model
% lambda == 1: Gaussian model
    res = zeros(size(x));
    for i = 1:amt
%         minindx = find(x > center + (i-1)*dist - 10*hw, 1, 'first');  % Does not decrease the runtime. Maybe in other circumstances 
%         maxindx = find(x < center + (i-1)*dist + 10*hw, 1, 'last');
        rel = BC(amt,i)/BC(amt, ceil((amt-1)/2)+1); % relative height of peak i 
%         res(minindx:maxindx) = res(minindx:maxindx) + gl(x(minindx:maxindx), hw, (center + (i-1)*dist), rel*height, lambda);
        res = res + gl(x, hw, (center + (i-1)*dist), rel*height, lambda);
    end
end

