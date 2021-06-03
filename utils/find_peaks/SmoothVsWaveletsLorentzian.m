% Script that compares segmented smoothing to wavelet denoise, for a
% signal consisting of 10 Lorentzian peaks whose widths vary randomly by
% a factor of 10 over the duration of the signal. The smooth widths (line
% 19) and the wavelet level (line 22) are chosen to optimize each method.
% Prints out the factor by which the wavelet method is slower.
% Requires SegmentedSmooth.m and wdenoise (from the Matlab wavelet
% toolbox). Compares slifing average, Savitsky-Golay, Fourier filter, 
% and two different wavelet types, including the coiflet2
% http://wavelets.pybytes.com/wavelet/coif2/ and
% Harr, http://wavelets.pybytes.com/wavelet/haar/
format short g
format compact
figure(1)
clf
% Initial settings (user settable)
increment=.5;
Noise=.5; % Amount of random noise added to the signal. (Change if desired) 
Baseline=2;
smoothwidths=11;
SmoothType=4;
ends=0;
Level=9; % Wavelet level

% load data file consisting of several peaks of increasing width
load DemoLorentzians.mat
demodata=[x' y']; % Assembles x and y vectors into data matrix
tic;ff=ifilter(x,y,0,.75,11,2,'Low-pass');fftoc=toc;
clf

subplot(2,3,1)
plot(x,y,'b')  % Graph the original signal in blue
thisaxis= [0 1200 -0.4 12];
axis(thisaxis)
xlabel('Original unsmoothed data')
title('SmoothVsWaveletLorentzian.m')
legend('Data with white noise')

subplot(2,3,2)
tic;sy=SegmentedSmooth(y,smoothwidths,SmoothType,ends);t1=toc;
NumSegments=length(smoothwidths);
plot(x,z,'b',x,sy,'r')  % Graph the smoothed signal in red
axis(thisaxis); % Use the same plotting axes as before.
xlabel(['RSD residual =  ' num2str(100.*std((z-sy)./max(z))) '%' ])
title(['Sliding average, width = ' num2str(smoothwidths) ])
legend('Noiseless data','Denoised noisy data')

subplot(2,3,3)
sgwidth=25;
tic;dny=ProcessSignal(x,y,0,sgwidth,4,0,0,0,0,0,0);t2=toc;
plot(x,z,'b',x,dny,'r')   % Graph the denoised signal 
axis(thisaxis);
title('Savitsky-Golay smooth, width=30')
title(['Savitsky-Golay smooth, width = ' num2str(sgwidth) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
SavitskyGolaySpeedPenalty=t2/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,4)
ffbandwidth=.9;
shape=6; % moderatley sharp cut-off rate
tic;ff=FouFilter(y,max(x),0,ffbandwidth,2,0);fftoc=toc;
plot(x,z,'b',x,ff,'r')  % Graph the Fourier filter signal 
axis(thisaxis);
title(['Fourier filter (low pass), bandwidth = ' num2str(ffbandwidth) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-ff)./max(z))) '%' ])
FourierFilterSpeedPenalty=fftoc/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,5)
tic;dny=wdenoise(y,'Wavelet','coif2');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title('Coiflets2 wavelet denoise')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
coif2WaveletSpeedPenalty=t2/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,6)
tic;dny=wdenoise(y,12,'Wavelet','haar');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')   % Graph the denoised signal 
axis(thisaxis);
title('Haar wavelet denoise')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
HaarWaveletSpeedPenalty=t2/t1
legend('Noiseless data','Denoised noisy data')