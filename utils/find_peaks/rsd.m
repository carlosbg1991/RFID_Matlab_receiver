function relstddev=rsd(x)
% Relative standard deviation of vector x
relstddev=std(x)./mean(x);