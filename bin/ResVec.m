%% Reduces dimension of optimizing routines, UNUSED
function R2 = ResVec(R, lenX, numvar)
    if length(R) == sum(lenX) % No constraints
        R2 = zeros(length(lenX)*numvar);
    else                      % Some constraints
        R2 = zeros(length(lenX)*numvar+length(R)-sum(lenX));
    end
    ind = 0;
    for i = 1:length(lenX)
        Rtemp = R(ind+1:ind+lenX(i)).^2;
        R2(i) = sqrt(sum(Rtemp));
        ind = ind + lenX(i);
    end
    if length(R) > sum(lenX)
        R2(length(lenX)+1:end) = R(ind+1:end);
    end
end
