% Searches for most negative values in a moving window to deal with noise
function [xmin, ymin] = findMinima3(d, ~, d2, bw, epsi, delta)
    % Initialising vectors
    d2xmin = 0; % Supposed maxima in d, due to minima in d2
    xmin   = 0; % Real maxima in d, thanks to hints from d2xmin
    ymin   = 0; % Minimal points in d2 with d > eps
    
    for i = (bw+1):(length(d2)-bw)
        if d2(i) < delta
            cl = (d2(i) < d2(i-bw:i-1));
            cr = (d2(i) < d2(i+1:i+bw));
            if all(cl)&&all(cr)  % Checks if all entries are true
               d2xmin = [d2xmin i]; %#ok<AGROW>
            end
        end
    end
    d2xmin = d2xmin(2:end);
    lend = length(d);
    for ind = d2xmin
       mux = max(d(ind-bw:ind+bw));
       if mux >= epsi 
           k = find(mux == d(ind-bw:ind+bw))+ind-bw-1;
           if k == ind-bw || k == ind+bw % Prevention of shoulder peaks, runs to local max in unsmoothed data
               k = ind; % First looks for local minima, if no local min is found, goes until the peak is found (even outside window)
               while (k>1 && d(k-1) > d(k)) 
                   k = k-1;
               end
               while (k<lend && d(k+1)>d(k))
                   k = k+1;
               end
               while (k>1 && d(k-1) > d(k)) % Tests first condition and ignores the second one if k == 0 
                   k = k-1;
               end  
           end
           xmin = [xmin k]; %#ok<AGROW>
           ymin = [ymin ind]; %#ok<AGROW>
       end
    end

    temp = xmin; % So we don't look up d(0) accidentally if we've already set the value at xmin(i+1) to 0
    currind = 2;
    for i=3:length(xmin)
        if xmin(i)-xmin(currind) < bw
            if d(xmin(currind)) >= d(xmin(i))
                temp(i) = 0;
            else
                temp(currind) = 0;
                currind = i;
            end
        else
            currind = i;
        end
    end

    xmin  = temp; % Retrieve from auxiliary vector
    xmin  = sort(unique(xmin)); % '0' always in front
    xmin  = xmin(2:end); 
    ymin  = ymin(2:end);
    postcorrect = d(xmin); % There were still edgecases where this part became necessary
    biginds = (postcorrect >= epsi);
    xmin = xmin(biginds);
    % ymin = ymin(biginds); % Leave line out so you can still see old peaks in plots
end