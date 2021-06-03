function NoiseColorTest
% Noise Color Test, version 3 (January 2017), by toh@umd.edu
% The frequency distribution of noise, designated by noise color,
% substantially effects the ability of smoothing to reduce noise. The
% Matlab/Octave function “NoiseColorTest.m” compares the effect of a
% 20-point boxcar (unweighted sliding average) smooth on the standard
% deviation of white, pink, and blue noise, all of which have an original
% unsmoothed standard deviation of 1.0. Because smoothing is a low-pass
% filter process, it effects low frequency (pink and red) noise less, and
% effects high-frequency (blue and violet) noise more, than it does white
% noise. Figure 1 shows the five type of noises before smoothing, figure 3
% shown their frequency spectra, and figure 2 shows the five type of noises
% after smoothing. 
format compact
% User-variable constants:
x=1:1000;
Noise=1;
SmoothWidth=20;
SmoothType=1; % 1=rectangular (sliding average); 2=triangle; 3=Gaussian

WhiteNoiseArray=Noise*(randn(size(x))); % normal white constant noise
BlueNoiseArray=Noise*bluenoise(length(x)); % Blue noise
PinkNoiseArray=Noise*pinknoise(length(x)); % Pink (1/f) noise
RedNoiseArray=cumsum(WhiteNoiseArray)./(std(cumsum(WhiteNoiseArray)));
VioletNoiseArray=Noise*violetnoise(length(x)-1); % Violet noise
size(VioletNoiseArray)

SmoothedWhiteNoise=fastsmooth(WhiteNoiseArray,SmoothWidth,SmoothType);
SmoothedPinkNoise=fastsmooth(PinkNoiseArray,SmoothWidth,SmoothType);
SmoothedBlueNoise=fastsmooth(BlueNoiseArray,SmoothWidth,SmoothType);
SmoothedRedNoise=fastsmooth(RedNoiseArray,SmoothWidth,SmoothType);
SmoothedVioletNoise=fastsmooth(VioletNoiseArray,SmoothWidth,SmoothType);

disp('Standard Deviations:')
BeforeSmoothingSTD=stdev(WhiteNoiseArray)
SmoothedWhiteSTD=stdev(SmoothedWhiteNoise)
SmoothedPinkSTD=stdev(SmoothedPinkNoise)
SmoothedBlueSTD=stdev(SmoothedBlueNoise)
SmoothedRedSTD=stdev(SmoothedRedNoise)
SmoothedVioletSTD=stdev(SmoothedVioletNoise)

figure(1)
clf
subplot(2,3,2)
plot(x,WhiteNoiseArray)
axis([0 length(x) -3 3])
title('White Noise')
xlabel(['Standard deviation = ' num2str(stdev(WhiteNoiseArray))] )
subplot(2,3,3)
plot(x,PinkNoiseArray)
axis([0 length(x) -3 3])
title('Pink (1/f) Noise')
xlabel(['Standard deviation = ' num2str(stdev(WhiteNoiseArray))] )
subplot(2,3,4)
plot(x,BlueNoiseArray)
axis([0 length(x) -3 3])
title('Blue Noise')
xlabel(['Standard deviation = ' num2str(stdev(PinkNoiseArray))] )
subplot(2,3,5)
plot(x,RedNoiseArray)
title('Red Noise')
xlabel(['Standard deviation = ' num2str(stdev(RedNoiseArray))] )
axis([0 length(x) -3 3])
subplot(2,3,6)
plot(x,VioletNoiseArray)
title('Violet Noise')
xlabel(['Standard deviation = ' num2str(stdev(VioletNoiseArray))] )
axis([0 length(x) -3 3])

figure(2)
subplot(2,3,1)
plot(x,WhiteNoiseArray,'c')
hold on
switch SmoothType,
    case 1
        plot(x,rectanglepulse(x,length(x)/2,SmoothWidth),'r')
    case 2
        plot(x,triangle(x,length(x)/2,SmoothWidth),'r')
    case 3
        plot(x,gaussian(x,length(x)/2,SmoothWidth),'r')
end
axis([0 length(x) -3 3])
hold off
title('White noise (Red=smoothing funciton)')
xlabel(['Width = '  num2str(SmoothWidth)    '   Type = '  num2str(SmoothType) ] )
subplot(2,3,2)
plot(x,SmoothedWhiteNoise,'c')
axis([0 length(x) -3 3])
title('Smoothed White Noise')
xlabel(['Standard deviation = ' num2str(SmoothedWhiteSTD)] )
subplot(2,3,3)
plot(x,SmoothedPinkNoise,'c')
axis([0 length(x) -3 3])
title('Smoothed Pink (1/f) Noise')
xlabel(['Standard deviation = ' num2str(SmoothedPinkSTD)] )
subplot(2,3,4)
plot(x,SmoothedBlueNoise,'c')
axis([0 length(x) -3 3])
title('Smoothed Blue Noise')
xlabel(['Standard deviation = ' num2str(SmoothedBlueSTD)] )
subplot(2,3,5)
plot(x,SmoothedRedNoise,'c')
title('Smoothed Red Noise')
xlabel(['Standard deviation = ' num2str(SmoothedRedSTD)] )
axis([0 length(x) -3 3])
subplot(2,3,6)
plot(x,SmoothedVioletNoise,'c')
title('Smoothed Violet Noise')
xlabel(['Standard deviation = ' num2str(SmoothedVioletSTD)] )
axis([0 length(x) -3 3])

figure(3);
clf
subplot(2,3,2)
PlotFrequencySpectrum(x,WhiteNoiseArray',4,0,0);
title('WhiteNoiseArray')
xlabel('Log-log frequency spectrum')
subplot(2,3,3)
PlotFrequencySpectrum(x,PinkNoiseArray',4,0,0);
title('PinkNoiseArray')
xlabel('Log-log frequency spectrum')
subplot(2,3,4)
PlotFrequencySpectrum(x,BlueNoiseArray',4,0,0);
title('BlueNoiseArray')
xlabel('Log-log frequency spectrum')
subplot(2,3,5)
PlotFrequencySpectrum(x,RedNoiseArray',4,0,0);
title('RedNoiseArray')
xlabel('Log-log frequency spectrum')
subplot(2,3,6)
PlotFrequencySpectrum(x,VioletNoiseArray',4,0,0);
title('VioletNoiseArray')
xlabel('Log-log frequency spectrum')

function y=whitenoise(x)
% Random noise with white power spectrum with mean zero 
% and unit standard deviation, equal in length to x
% Tom O'Haver, 2008
% Example:
% model=gaussian([1:10deriv00],500,200);
% model=model+.1.*whitenoise(model);
% plot(model)
% 
y=randn(size(x));

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
ry=deriv1(y);
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
lft1=[1:(length(fy)/2)+1];
lft2=[(length(fy)/2):length(fy)];
ffilter1=ones(size(lft1))./(sqrt(lft1));
ffilter2=ones(size(lft2))./(sqrt(lft2));
ffilter=[ffilter1,ffilter2];
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter(1:length(fy));  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Inverse transform to recover filtered signal 'ry'
ry=((ry-mean(ry))./stdev(ry)); % Normalize to unit standard deviation

function vy = violetnoise(x)
% Random noise with high-frequency weighted power spectrum with mean zero
% and unit standard deviation. Length equal to the length of x
y=0:x;
for n=1:2:length(y)-1,
    rn=abs(randn());
    y(n)=rn;
    y(n+1)=-rn;
end
vy=(y-mean(y))./stdev(y); % Normalize to unit standard deviation

function d=deriv1(a)
% First derivative of vector using simple adjacent differences.
% Example: deriv([1 1 1 2 3 4]) yeilds [0 0 0 1 1 1]
%  T. C. O'Haver, 2013.
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
for j = 2:n;
  d(j)=a(j)-a(j-1);
end
