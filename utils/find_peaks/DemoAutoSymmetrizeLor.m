% Demonstration of symmetrization of an exponentially modified
% Lorentzian (EML) by weighted addition of the first derivative
% using the AutoSymmetrize function. Requires AutoSymmetrize and
% lorentzian functions.
format compact
format short g
clear
A=1; % peak height
s=100; % half-width of original Lorentzian, in time units
tau=100; % Exponential broadening factor, in time units
mu=900; % peak center, in time units
Noise=.001; % random white noise
SmoothWidth=7; % Smooth width to reduce effect of noise
plots=0;  % 1 = show plots within AutoSymmetrize, 0 = no plots within AutoSymmetrize
deltat=5; % sampling interval, in time units
t=0:deltat:2500; % time axis

% Create an exponentially broadened Gaussian
EML=A.*lorentzian(t,mu,s); %  original Lorentzian
taupoints=tau/deltat; % tau expressed in number of points
EMLb=ExpBroaden(EML',-taupoints)'; %  broadened Lorentzian
EMLbn=EMLb+Noise.*randn(size(t)); %  noisy broadened Lorentzian

% Call the AutoSymmetrize function
[EMLSymm,EstimatedTau]=AutoSymmetrize(t,EMLbn,SmoothWidth,plots);

% Optionally plot the results
if plots
else
    figure(1)
    clf
    plot(t,EMLbn,t,EMLSymm,t,max(EMLSymm).*EML,'k.')
    xlabel('Black dots: Original Lorentzian   Blue: Observed exponentially broadened peak.   Red: Symmetrized peak')
    title(['Automated determination of exponential broadening factor, tau = '  num2str(EstimatedTau) '.'])
end
disp(' ')
MeasuredTau=EstimatedTau
OriginalArea=trapz(t,EML) % Trapeziodal area of original peak
SymmetrizedArea=trapz(t,EMLSymm) % Trapeziodal area of symmetrized peak