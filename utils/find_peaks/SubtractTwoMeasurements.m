% Demonstration of measuring the noise and signal-to-noise ratio of a
% stable waveform, by subtracting two measurements. The signal must be
% stable between measurements. See
% https://terpconnect.umd.edu/~toh/spectrum/SignalsAndNoise.html or
% https://terpconnect.umd.edu/~toh/spectrum/CaseStudies.html#SNR
% Tom O'Haver, December 2016
x=2:.001:7;
Noise=1;
% n1=Noise.*randn(size(x));  % White noise in measurement 1
n1=Noise.*pinknoise(length(x)); % Pink noise in measurement 1
y1=100.*exp(-(x-2).^2)+100.*exp(-(x-4).^2)+30.*exp(-(x-6).^2);% Pure signal 1
m1=y1+n1; % Measurement 1, with added noise
% n2=Noise.*randn(size(x)); % White noise in measurement 2
n2=Noise.*pinknoise(length(x)); % Pink noise in measurement 2
y2=100.*exp(-(x-2).^2)+100.*exp(-(x-4).^2)+30.*exp(-(x-6).^2);% Pure signal 2
m2=y2+n2;% Measurement 2, with independent added noise (not the same as n1)
subplot(3,1,1);plot(x,m1)
title('First measurement: m1')
subplot(3,1,2);plot(x,m2)
title('Second measurement: m2')
subplot(3,1,3);plot(x,m1-m2)
title('Difference between m1 and m2')
MeasuredNoise=sqrt((stdev(m1-m2).^2)/2); 
SNR=max(y1)./MeasuredNoise; % SNR = signal-to-noise ratio
disp(' ')
disp(['Actual noise standard deviation was ' num2str(std([n1 n2]))])
disp(['Measured noise in m1-m2: ' num2str(std(m1-m2)) ] )
disp(['Calculated noise in signal: ' num2str(MeasuredNoise) ] )
disp(['Measured signal-to-noise ratio is ' num2str(SNR)])

% Derivation
%  signal=noiseless underlying signal
%  n1=random noise present in measurement 1, average value = zero
%  n2=random noise present in measurement 2 (Separate sample, same population)
%  sdn=standard deviation of n1 or n2 (same population, on average equal)
%  measurement 1: m1=signal+n1   assuming noise is added to signal
%  measurement 2: m2=signal+n2 (signal is exactly the same in both)
%  difference=m1-m2, so that signal cancels out exactly, leaving noise
% 
% (See first line of https://terpconnect.umd.edu/~toh/spectrum/ErrorPropagation.pdf)
%  Rules for error propagation say that when two random variables are added 
%  or substrated, the standard deviation of the result is the square root 
%  of the sum of the squares of the individual standard deviations:
% 
%  std(m1-m2) = sqrt(std(sdn)^2 + std(sdn)^2)
% 
% Square both sides:
%  std(m1-m2)^2 = std(sdn)^2 + std(sdn)^2
% 
% Combine terms on right side:
%  std(m1-m2)^2 = 2*(std(sdn)^2)
% 
% Divide both sides by 2:
%  (std(m1-m2)^2)/2 = std(sdn)^2)
% 
% Take square root of both sides
%  sqrt((std(m1-m2)^2)/2) = std(sdn)
