% Overhead function for 'mysg.m'
function[y,y_abl1,y_abl2,y_abl3]= mysgfilt(d, N, x)
%% Preparation phase
    L      = length(x);
    y      = zeros(1,L);
    y_abl1 = zeros(1,L);
    y_abl2 = zeros(1,L);
    y_abl3 = zeros(1,L);
    if mod(N,2) == 0
        warning('N must be odd! Returned zeros.')
        return
    end
    M = (N-1)/2;
    if d >= N
        warning('d must be strictly smaller than N! For proper smoothing results, choose d to be greatly smaller than N. Returned zeros.')
        return
    end

    [H,~] = mysg(d,N); % H contains convoultion vectors
%% Calculate first M values
    polco = zeros(d+1,1);
    firx = x(1:2*M+1);
    for i = 0:d
        polco(i+1) = H(d+1-i,:)*firx;
    end
    for i=1:M
        y(i) = polyval(polco, -M+i-1);
    end
    if d >= 1
        polco = polyder(polco);
        for i=1:M
            y_abl1(i) = polyval(polco, -M+i-1);
        end
        if d >= 2
            polco = polyder(polco);
            for i=1:M
                y_abl2(i) = polyval(polco, -M+i-1);
            end
            if d >= 3
                polco = polyder(polco);
                for i=1:M
                y_abl3(i) = polyval(polco, -M+i-1);
                end
            end
        end
    end
%% Calculate main values
    X = zeros(N,L-2*M);
    for j=1:size(X,2)
        X(:,j) = x(j:j+2*M);
    end
    y(M+1:L-M) = H(1,:)*X;
    if d >= 1
        y_abl1(M+1:L-M) = 1*H(2,:)*X;
        if d >= 2
            y_abl2(M+1:L-M) = 2*H(3,:)*X;
            if d >= 3
                y_abl3(M+1:L-M) = 6*H(4,:)*X;
            end
        end
    end
%% Calculate last M values
    polco = zeros(d+1,1);
    lasx  = x(L-2*M:L);
    for i = 0:d
        polco(i+1) = H(d+1-i,:)*lasx;
    end
    for i=1:M
        y(L-M+i) = polyval(polco, i);
    end
    if d >= 1
        polco = polyder(polco);
        for i=1:M
            y_abl1(L-M+i) = polyval(polco, i);
        end
        if d >= 2
            polco = polyder(polco);
            for i=1:M
                y_abl2(L-M+i) = polyval(polco, i);
            end
            if d >= 3
                polco = polyder(polco);
                for i=1:M
                y_abl3(L-M+i) = polyval(polco, i);
                end
            end
        end
    end
    y      = y';
    y_abl1 = y_abl1';
    y_abl2 = y_abl2';
    y_abl3 = y_abl3';
end