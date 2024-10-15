%% Function doing the actual Matrix-calculations necessary for SG
function [H, A] = mysg(d, N) 
    M = (N-1)/2;
    A = zeros(N,d+1);
    for m=-M:M
       for i=0:d
          A(m+M+1, i+1) = m^i;
       end
    end
    
    F = A' * A;
    H = (inv(F))*A'; % 'H(i,:)*x(x0-M:x0+M)' = 'i-1'-th derivative of the SG-Polynomial at x = x0.
end