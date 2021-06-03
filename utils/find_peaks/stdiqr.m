function s=stdiqr(x)
% function s=stdiqr(x) returns an estimate of the standard deviation of
% vector x based on the inter-quartile range (irq), which is more robust to
% large outliers that the std function.
% Examples:
% std(randn(1,1000000)) returns a number very close to 1.000
% stdiqr(randn(1,1000000)) returns a number very close to 1.000
% std([2000 randn(1,1000000)]) returns 2.23 because of a single outlier
% stdiqr([2000 randn(1,1000000)]) returns number close to 1.
sizex=size(x);
if sizex(1)<sizex(2),x=x';end
size(x)
s=IQrange(x);

function b=IQrange(a)
% b = IQrange(a) returns the interquartile range divided by 1.34896, of the
% values in a, which is an estimate of the standard deviation of a without
% ouliers. If a is a vector, b is the difference between the 75th and 25th
% percentiles of a, divided by 1.34896. If b is a matrix, b is a row vector
% containing the interquartile range of each column of a, divided by
% 1.34896. The factor 1.34896 makes this quantity equal on average to the
% standard deviation (SD) of normal distributions, and thus IQrange is a
% better estimate of the standard deviation without ouliers for a normal
% distribution.
%  T. C. O'Haver, 2017
% Example:
%  a=randn(10000,5);
%  IQrange(a)
% Divide by 1.34896 to get an estimate of the standard deviation withput
% outliers
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