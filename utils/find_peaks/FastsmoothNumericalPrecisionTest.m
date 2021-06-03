% Demonstration of numerical precision of fastsmooth.m Typical result for
% 4000-point Gaussian smooth applied to a 100,000 point signal:
% RSDdiffCenter =
%       3.59274620923292e-17
% RSDdiffAll =
%       2.69007007479015e-06

data=1:100000;
SmoothWidth=4000;
SmoothType=3;
sdata=fastsmooth(data,SmoothWidth,SmoothType,1);
diff=data-sdata;
clf
plot(diff);
% Compute relative standard deviation (RSD) of the center portion of the
% signal and of the entire signal includinf the ends.
RSDdiffCenter=std(diff(SmoothWidth*SmoothType:length(data)-SmoothWidth*SmoothType))./max(data)
RSDdiffAll=std(diff(1:length(data)))./max(data)