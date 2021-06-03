% Script that compares segmented smoothing to wavelet denoise, for a
% signal consisting of peaks whose widths increase gradually by
% a factor of 10 over the duration of the signal. The smooth widths (line
% 16) and the wavelet level (line 19) are chosen to optimize each method.
% Prints out the factor by which the wavelet method is slower.
% Requires SegmentedSmooth.m and wdenoise (from the Matlab wavelet
% toolbox). Compares 4 different wavelet types.
format short g
format compact
figure(1)
clf
% Initial settings (user settable)
increment=.5;
Noise=.5; % Amount of random noise added to the signal. (Change if desired) 
Baseline=2;
smoothwidths=411;
SmoothType=1;
ends=0;
Level=9; % Wavelet level

% load data file consisting of several peaks of increasing width
load DemoSquareWave.mat

demodata=[x' y']; % Assembles x and y vectors into data matrix
tic;ff=ifilter(x,y,0.03,0.01,1,2,'Low-pass');fftoc=toc;
clf
subplot(2,3,1)
plot(x,y,'b')  % Graph the original signal in blue
thisaxis=axis;
xlabel('Original unsmoothed data')
title('SmoothVsWaveletsSquareWave.m')
legend('Data with white noise')

subplot(2,3,2)
tic;sy=SegmentedSmooth(y,smoothwidths,SmoothType,ends);t1=toc;
NumSegments=length(smoothwidths);
plot(x,z,'b',x,sy,'r')  % Graph the smoothed signal in red
axis(thisaxis); % Use the same plotting axes as before.
xlabel(['RSD residual =  ' num2str(100.*std((z-sy)./max(z))) '%' ])
title(['Sliding average, width =  ' num2str(smoothwidths) ])
legend('Noiseless data','Denoised noisy data')

subplot(2,3,3)
tic;dny=ProcessSignal(x,y,0,300,4,0,0,0,0,0,0);t2=toc;
plot(x,z,'b',x,dny,'r')   % Graph the denoised signal 
axis(thisaxis);
title('Savitsky-Golay smooth, width=300')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
SavitskyGolaySpeedPenalty=t2/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,4)
plot(x,z,'b',x,ff,'r')   % Graph the Fourier filter signal 
axis(thisaxis);
title('Fourier filter (low pass)')
xlabel(['RSD residual =  ' num2str(100.*std((z-ff)./max(z))) '%' ])
FourierFilterSpeedPenalty=fftoc/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,5)
tic;dny=wdenoise(y,'Wavelet','coif2');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')   % Graph the denoised signal 
axis(thisaxis);
title('coif2 wavelet denoise')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
coif2WaveletSpeedPenalty=t2/t1 % How much slower this function is
legend('Noiseless data','Denoised noisy data')

subplot(2,3,6)
tic;dny=wdenoise(y,15,'Wavelet','haar');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')   % Graph the denoised signal 
axis(thisaxis);
title('haar wavelet denoise')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
HaarWaveletSpeedPenalty=t2/t1
legend('Noiseless data','Denoised noisy data')