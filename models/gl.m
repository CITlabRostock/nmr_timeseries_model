function [res] = gl(x,a,b,c,lambda)
if lambda==0 % Lorentz
    res=(c*a^2)./((x-b).^2+a^2);
elseif lambda==1 % Gauss
    res=c.*exp((-log(2)/(a^2)).*(x-b).^2);
else %gauss_lorentz
    y = (x-b).^2;               % 2 Vector_ops
    res=(lambda*c).*exp((-log(2)/(a^2)).*y) + ((1-lambda)*c*a^2)./(y+a^2); % 3 Vec_ops for first summand, 2 for second, 1 for adding
end % 8 Vector operations in total

