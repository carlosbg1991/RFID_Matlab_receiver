function y=propnoise(x)
% Random noise whose amplitude is proportional to the amplitude of signal
% zwith white power spectrum with mean zero and unit standard deviation,
% equal in length to x.
% Tom O'Haver, 2012
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*propnoise(model);
% plot(model)
% 
z=x.*randn(size(x));
sz=std(z);
y=z./sz;
