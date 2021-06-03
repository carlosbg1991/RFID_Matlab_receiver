function by=bluenoise(n)
% Random noise with blue power spectrum with mean zero 
% and unit standard deviation. n is number of points.
% Tom O'Haver, 2008
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*bluenoise(length(model));
% plot(model)
x=1:n;
y=randn(size(x));  % Random normally-distributed white noise
ry=deriv(y);
by=(ry./stdev(ry)); % Normalize to unit standard deviation

function stddev=stdev(a)
% Octave and Matlab compatible standard deviation function
sa=size(a);
reshape(a,1,length(a));
if sa(1)>sa(2),
  stddev=std(a);
else
  stddev=(std(a'));
end;