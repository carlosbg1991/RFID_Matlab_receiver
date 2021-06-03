% Demonstration of the effect of amplitude noise, frequency noise, and
% phase noise on the ensemble averaging of a sine waveform signal. Shows
% that (a) ensemble averaging reduces the white noise in the signal but not
% the frequency or phase noise; (b) ensemble averaging the Fourier
% transform has the same effect as ensemble averaging the signal itself;
% and (c) the effect of phase noise is reduced if the power spectra are
% ensemble averaged. 
increment=.1; % difference between adjacent x values
x=0:increment:10*pi;
Amplitude=1;
WhiteNoise=1; % Point-to-point changes in the amplitude of each signal
FreqNoise=.1; % Signal-to-signal changes in the frequency of the sine wave
PhaseNoise=1; % Signal-to-signal changes in the phase of the sine wave
NumSignals=200; % Number of signal averaged
sumy=0;
sumfft=0;
sumpower=0;
clear signalmatrix
sy=Amplitude.*sin(x);
for n=1:NumSignals,
    y=Amplitude.*sin(x.*(1+randn()*FreqNoise)+randn().*PhaseNoise);
    noisyy=y+WhiteNoise.*randn(size(x));
    sumy=sumy+noisyy;
    sumfft=sumfft+fft(noisyy);
    sumpower=sumpower+real(sqrt(fft(noisyy) .* conj(fft(noisyy))));
    signalmatrix(n,:)=noisyy;
end
Averagey=sumy./NumSignals;
Averagefft=sumfft./NumSignals;
AveragePower=sumpower./NumSignals;
figure(1)
subplot(3,1,1)
plot(x,Averagey,x,sy,'-.')
title('Ensemble averaged signals (Dotted line is original noiseless signal)')
subplot(3,1,2)
plot(x,ifft(Averagefft),x,sy,'-.')
title('Inverse Fourier transform of ensemble averaged Fourier transforms')
subplot(3,1,3)
plot(x,AveragePower)
plot(x,ifft(AveragePower),x,sy,'-.')
AXIS([min(x) max(x) min(sy) max(sy)])
title('Inverse Fourier transform of ensemble averaged power spectra')
PowerSum=sum(ifft(AveragePower));
PowerMean=mean(ifft(AveragePower));