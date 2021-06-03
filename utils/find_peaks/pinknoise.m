function ry=pinknoise(n)
% Random noise with pink (1/f) power spectrum with mean zero 
% and unit standard deviation. n is number of points.
% Tom O'Haver, 2008
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*pinknoise(length(model));
% plot(model)
x=[1:n];
y=randn(size(x));  % Random normally-distributed white noise
% Fourier filter 
fy=fft(y); % Compute Fourier transform of signal y
% Compute filter shape
lft1=[1:(length(fy)/2)-1];
lft2=[(length(fy)/2):length(fy)];
ffilter1=ones(size(lft1))./(sqrt(lft1));
ffilter2=ones(size(lft2))./(sqrt(lft2));
ffilter=[ffilter1,ffilter2];
if length(fy)>length(ffilter),ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter(1:length(fy));  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Inverse transform to recover filtered signal 'ry'
% figure(2);plot(real(ry))
ry=((ry-mean(ry))./stdev(ry)); % Normalize to zero mean and unit standard deviation

function stddev=stdev(a)
% Octave and Matlab compatible standard deviation function
sa=size(a);
reshape(a,1,length(a));
if sa(1)>sa(2),
  stddev=std(a);
else
  stddev=(std(a'));
end;