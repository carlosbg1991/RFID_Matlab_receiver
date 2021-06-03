function b=IQrange(a)
% b = IQrange(a) returns the interquartile range of the values in "a",
% divided by 1.34896, which is an estimate of the standard deviation of "a"
% without ouliers. If "a" is a vector, b is the difference between the 75th
% and 25th percentiles of a, divided by 1.34896. If "a" is a matrix, b is a
% row vector containing the interquartile range of each column of a,
% divided by 1.34896. The factor 1.34896 makes this quantity equal on
% average to the standard deviation (SD) of normal distributions, and thus
% IQrange is a better estimate of the standard deviation without ouliers
% for a normal distribution.
%  T. C. O'Haver, 2017
% Example:
%  a=randn(1000,1); % Create randome numbers with standard deviation of 1.
%  std(a) % Display standard deviation
%  a(500)=100; % Add single outlier
%  std(a) % Shows that standard deviation is larger with single outlier.
%  IQrange(a)  % Shows that IQrange is not effected by outlier.

mina=min(a);
sizea=size(a);
NumCols=sizea(2);
b=zeros(size(a));
for n=1:NumCols
    b(:,n)=a(:,n)-mina(n);
end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));
b=b./1.34896;