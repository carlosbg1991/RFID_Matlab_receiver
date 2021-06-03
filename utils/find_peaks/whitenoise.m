function y=whitenoise(x);
% Random noise with white power spectrum with mean zero 
% and unit standard deviation, equal in length to x
% Tom O'Haver, 2008
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*whitenoise(model);
% plot(model)
% 
y=randn(size(x));
