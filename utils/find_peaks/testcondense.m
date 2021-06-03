% Demonstration of the effect of boxcar averaging using the consense.m
% function to reduce noise without changing the noise color. Shows that the
% boxcar reduces the measured noise, removing the high frequency
% components but has little effect on the the peak parameters.
% Least-squares curve fitting on the condensed data is faster and results
% in a lower fitting error, but no more accurate measurement of peak parameters. 
cf=25; % Condense factor = number of points in boxcar average (try 2-100)
noise=0.1; % Percent random white noise
x=1:.02:100;
S=lorentzian(x,45,10) + lorentzian(x,60,10); % signal vector
N=noise.*randn(size(x)); % noise vector
y=S+N;
disp(' ')
OriginalSignalToNoiseRatio=max(S)./stdev(N);
disp(['Original Signal-To-Noise Ratio =' num2str(OriginalSignalToNoiseRatio) ] )
cx=condense(x,cf); % cx = condensed x vector
cy=condense(y,cf); % cy = condensed y vector
cn=condense(N,cf); % cy = condensed noise vector
CondensedSignalToNoiseRatio=max(cy)./stdev(cn);
disp(['Condensed Signal-To-Noise Ratio =' num2str(CondensedSignalToNoiseRatio) ] )

figure(1)
clf
plot(x,y,'.c',cx,cy,'.r')
title('Comparison of original signal (cyan) and condensed signal (red) ')
xlabel('x')
ylabel('y (cyan) and condensed y (red)')

figure(2)
clf
PS=PlotFrequencySpectrum(x', y', 4, 0, 0);
PSc=PlotFrequencySpectrum(cx', cy', 4, 0, 0);
% Adjust the amplitudes of the frequency spectra to compensate for
% different lengths of original and condensed signal
loglog(PS(:,1),(PS(:,2)./(length(x).^2)),'c');
hold on;
loglog(PSc(:,1),(PSc(:,2)./(length(cx).^2)),'r')
hold off
ylabel('Amplitude')
xlabel('Frequency')
title('Comparison of frequency spectra of original (cyan) and condensed (red) signals')

figure(3)
disp('2-peak Voigt model fit to original signal')
tic
[PeakTable,GoodnessOfFit]=peakfit([x;y],50,62,2,20,2)
toc

disp(' ')
disp('2-peak Voigt model fit to condensed signal')
figure(4)
tic
[PeakTable,GoodnessOfFit]=peakfit([cx;cy],50,62,2,20,2)
toc