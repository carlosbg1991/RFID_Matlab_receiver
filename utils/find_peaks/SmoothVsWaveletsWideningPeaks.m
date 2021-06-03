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
Level=9; % Wavelet level

demodata=[x' y']; % Assembles x and y vectors into data matrix
clf
subplot(2,3,1)
plot(x,y,'b')  % Graph the original signal in blue
thisaxis=axis;
xlabel('Original unsmoothed data')
title('SmoothVsWaveletsWideningPeaks.m')
legend('Data with white noise')

subplot(2,3,2)
tic;sy=SegmentedSmooth(y,smoothwidths,SmoothType,ends);t1=toc;
NumSegments=length(smoothwidths);
plot(x,z,'b',x,sy,'r')  % Graph the smoothed signal in red
axis(thisaxis); % Use the same plotting axes as before.
xlabel(['RSD residual =  ' num2str(100.*std((z-sy)./max(z))) '%' ])
title(['Segmented smooth widths =  ' num2str(smoothwidths) ])
legend('Noiseless data','Denoised noisy data')

subplot(2,3,3)
sgwidth=80;
tic;dny=ProcessSignal(x,y,0,sgwidth,4,0,0,0,0,0,0);t2=toc;
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title(['Savitsky-Golay smooth width = ' num2str(sgwidth) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
SavitskyGolaySpeedPenalty=t2/t1 % How much slower this function is

subplot(2,3,4)
shape=3;
bandwidths=[20 19 17 15 10];
cf=ones(1,length(bandwidths));
tic;sff=SegmentedFouFilter(y,2*pi,cf,bandwidths,shape,0);t2=toc;
plot(x,z,'b',x,sff,'r')  % Graph the Fourier filter signal 
axis(thisaxis);
title(['Segmented Fourier filter, bandwidths =  ' num2str(bandwidths) ])
xlabel(['RSD residual =  ' num2str(100.*std((z-sff)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
FourierFilterSpeedPenalty=t2/t1 % How much slower this function is

subplot(2,3,5)
tic;dny=wdenoise(y,'Wavelet','coif2');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title('coif2 wavelet denoise')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
coif2WaveletSpeedPenalty=t2/t1 % How much slower this function is

subplot(2,3,6)
tic;dny=wdenoise(y,Level,'Wavelet','haar');t2=toc;plot(x,dny,'r')
plot(x,z,'b',x,dny,'r')  % Graph the denoised signal 
axis(thisaxis);
title('haar wavelet denoise')
xlabel(['RSD residual =  ' num2str(100.*std((z-dny)./max(z))) '%' ])
legend('Noiseless data','Denoised noisy data')
HaarWaveletSpeedPenalty=t2/t1