% Script that compares segmented smoothing to wavelet denoise, for a
% signal consisting of several peaks whose widths increase gradually by
% a factor of 10 over the duration of the signal. The smooth widths (line
% 16) and the wavelet level (line 19) are chosen to optimize each method.
% Prints out the factor by which the wavelet method is slower.
% Requires SegmentedSmooth.m and wdenoise (from the Matlab wavelet
% toolbox). Compares 4 different wavelet types.
format short g
format compact
figure(1)
clf
% load data file consisting of several peaks of increasing width
load WideningPeaks.mat % load one data set created by DemoSegmentedSmooth.m

% Initial settings (user settable)
smoothwidths=[25 35 55 60];
SmoothType=3;
ends=0;
Level=floor(log2(length(x))); % Wavelet level

demodata=[x' y']; % Assembles x and y vectors into data matrix
clf
subplot(2,3,1)
plot(x,y,'b')  % Graph the original signal in blue
thisaxis=axis;
xlabel('Original unsmoothed data')
title('WaveletsComparisonWideningPeaks.m')
legend('Data with white noise')

subplot(2,3,2)
tic;dny=wdenoise(y,Level,'DenoisingMethod','BlockJS');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title(['BlockJS wavelet denoise, level = ' num2str(Level) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
BlockJSWaveletSpeedPenalty=t2/t1 % How much slower this function is

subplot(2,3,3)
tic;dny=wdenoise(y,Level,'Wavelet','bior5.5');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title(['bior5.5 wavelet denoise, level = ' num2str(Level) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
bior5dot5SpeedPenalty=t2/t1 % How much slower this function is

subplot(2,3,4)
tic;dny=wdenoise(y,Level,'Wavelet','coif2');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title(['coif2 wavelet denoise, level = ' num2str(Level) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
coif2WaveletSpeedPenalty=t2/t1 % How much slower this function is

subplot(2,3,5)
tic;dny=wdenoise(y,Level,'Wavelet','sym8');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title(['sym8 wavelet denoise, level = ' num2str(Level) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
sym8WaveletSpeedPenalty=t2/t1 % How much slower this function is

subplot(2,3,6)
tic;dny=wdenoise(y,Level,'Wavelet','db4');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title(['db4 wavelet denoise, level = ' num2str(Level) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
db4WaveletSpeedPenalty=t2/t1