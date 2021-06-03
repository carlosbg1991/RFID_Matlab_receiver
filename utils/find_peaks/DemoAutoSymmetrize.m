% Demonstration of symmetrization of an exponentially modified
% Gaussian (EMG) by weighted addition of the first derivative
% using the AutoSymmetrize function. Requires AutoSymmetrize and
% expgaussian2 functions.
format compact
format short g
clear
% User-modifiable parameters
A=1; % peak height
s=100; % standard deviation of Gaussian peak
tau=100; % Exponential broadening factor
mu=1000; % peak center
Noise=.01; % random white noise
SmoothWidth=3; % Smooth width to reduce effect of noise
SmoothMode=3;
plots=0;  % 1 = show plots, 0 = no plots
deltat=10; % sampling interval
t=0:deltat:2500; % time axis

% Create an exponentially broadened Gaussian
EMG=A.*expgaussian2(t,mu,s,tau)+Noise.*randn(size(t));

% Call the AutoSymmetrize function
[EMGSymm,EstimatedTau]=AutoSymmetrize(t,EMG,SmoothWidth,plots);

% Print out tau and compare peak areas
MeasuredTau=EstimatedTau
OriginalArea=trapz(t,EMG) % Trapeziodal area of original peak
SymmetrizedArea=trapz(t,EMGSymm) % Trapeziodal area of symmetrized peak
SmoothedEMG=fastsmooth(EMG,SmoothWidth,SmoothMode,1);
disp(' ')
if plots
else
    figure(1)
    clf
    plot(t,SmoothedEMG,t,EMGSymm,t,max(EMGSymm).*gaussian(t,mu,2*sqrt(2*log(2))*s),'k.')
    xlabel('Black dots: Original Gaussian   Blue: Observed exponentially broadened peak.   Red: Symmetrized peak')
    title(['Automated determination of exponential broadening factor, tau = '  num2str(EstimatedTau) '.'])
end

MeasuredTau=EstimatedTau
OriginalArea=trapz(t,EMG) % Trapeziodal area of original peak
SymmetrizedArea=trapz(t,EMGSymm) % Trapeziodal area of symmetrized peak