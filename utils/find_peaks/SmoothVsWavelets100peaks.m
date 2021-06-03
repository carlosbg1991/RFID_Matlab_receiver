% Script that compares segmented smoothing to wavelet denoise, for a ecg
% signal. Methods are compared by computing the STD of the difference
% between the noiseless signal and the denoised noisy signal. Prints out
% the factors by which the various methods are slower that a sliding
% average. Compares sliging average, Savitsky-Golay, Fourier filter, and
% two different wavelet types, including the coiflet2
% http://wavelets.pybytes.com/wavelet/coif2/ and
% the Harr, http://wavelets.pybytes.com/wavelet/haar/
% Requires fastsmooth.m, ifilter.m, ProcessSignal,m, wdenoise
% (from the Matlab wavelet toolbox), and ecg.mat. Compares sliding average,

format short g
format compact
figure(1)
clf
% Initial settings (user settable)
increment=.5;
Noise=.3; % Amount of random noise added to the signal. (Change if desired) 
Baseline=2;
smoothwidths=5;
SmoothType=3;
ends=0;
Level=9; % Wavelet level

% load data file consisting of several peaks of increasing width
load HundredPeaks % A "clean" real ECG signal with little noise.

clf

subplot(2,3,1)
plot(x,y,'b')  % Graph the original signal in blue
thisaxis=axis;
xlabel('Original unsmoothed data')
title('SmoothVsWavelets100peaks.m')
legend('Data with white noise')

subplot(2,3,2)
tic;sy=fastsmooth(y,smoothwidths,SmoothType,ends);t1=toc;
NumSegments=length(smoothwidths);
plot(x,z,'b',x,sy,'r')  % Graph the smoothed signal in red
axis(thisaxis); % Use the same plotting axes as before.
xlabel(['RSD residual =  ' num2str(100.*std((z-sy)./max(z))) '%' ])
title(['Sliding average, width = ' num2str(smoothwidths) ])
legend('Noiseless data','Denoised noisy data')

subplot(2,3,3)
sgwidth=9;
tic;dny=ProcessSignal(x,y,0,sgwidth,4,0,0,0,0,0,0);t2=toc;
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title(['Savitsky-Golay smooth, width = ' num2str(sgwidth) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
SavitskyGolaySpeedPenalty=t2/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,4)
ffbandwidth=.3;
shape=6; % moderatley sharp cut-off rate
tic;ff=FouFilter(y,max(x),0,ffbandwidth,2,0);fftoc=toc;
plot(x,z,'b',x,ff,'r')  % Graph the Fourier filter signal 
axis(thisaxis);
title('')
title(['Fourier filter (low pass), bandpass = ' num2str(ffbandwidth) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-ff)./max(z))) '%' ])
FourierFilterSpeedPenalty=fftoc/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,5)
tic;dny=wdenoise(y,'Wavelet','coif2');t2=toc;
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title('Coiflets2 wavelet denoise')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
coif2WaveletSpeedPenalty=t2/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,6)
tic;dny=wdenoise(y,12,'Wavelet','haar');t2=toc;
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title('Haar wavelet denoise')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
HaarWaveletSpeedPenalty=t2/t1
legend('Noiseless data','Denoised noisy data')