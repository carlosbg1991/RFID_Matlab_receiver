function y=sqrtnoise(x)
% Random noise whose amplitude is proportional to the square root of the
% amplitude of signal with white power spectrum with mean zero and unit
% standard deviation, equal in length to x.
% Tom O'Haver, 2012
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*sqrtnoise(model);
% plot(model)
% 
z=sqrt(abs(x)).*randn(size(x));
sz=std(z);
y=z./sz;
