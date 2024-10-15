% Searches for most negative values in a moving window to deal with noise,
% UNUSED Predecessor of findMinima3.m
function [xmin, ymin] = findMinima2(d, ~, d2, bw, epsi, delta)
    %Initialising vectors
    xmin = 0;
    
    for i = (bw+1):(length(d2)-bw)
        if d2(i) <= delta
            cl = (d2(i) < d2(i-bw:i-1));
            cr = (d2(i) < d2(i+1:i+bw));
            if all(cr)&&all(cl)  % Checks if all entries are true
               mux = max(d(i-bw:i+bw));   % Might result in shoulderpeaks
               if mux >= epsi
                   k = find(mux == d(i-bw:i+bw)); % Index in small array
                   k = k+i-bw-1;                  % Index in large array
                   xmin = [xmin k]; %#ok<AGROW>
               end
            end
        end
    end
    xmin = xmin(2:end);  % Delete initial point from vectors
    xmin = unique(xmin); % Because of this badly thought out method this might become necessary
    ymin = d2(xmin);
end
