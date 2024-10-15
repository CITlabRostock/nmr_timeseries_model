function wstd = wasserstein(x, df1, df2)
% Calculates the Wasserstein-metric of two probability density functions.
% (If x is the approximately equidistant support of the two pdfs. Can deal
% with noise and even negative values
if ~(all(size(df1) == size(df2)) && all(size(x) == size(df1)))
    error('Dimensions of pdfs are unequal. Cannot calculate Wasserstein-metric')
end
% Norming
df1 = df1./sum(df1);
df2 = df2./sum(df2);
s1 = 0;
s2 = 0;
wstd = 0;
md = mean(diff(x)); % This suffices, since x is equidistant
for i=1:length(df1)
    s1 = s1 + df1(i); % Cumulative density function at point 'i'
    s2 = s2 + df2(i);
    wstd = wstd + abs(s1-s2)*md; % Calculates integral('abs(cdf1-cdf2))') in discrete setting in linear time.
end