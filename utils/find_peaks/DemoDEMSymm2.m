% Demonstration of symmetrization of a double exponentially modified
% Gaussian (EMG) by two successive applications of weighted addition of the
% first derivative, with weighting factors equal to the two taus. Shows
% that the peak width is reduced, the peak made more symmetrical, and the
% area is preserved. This demo scripts creates two overlapping double
% exponential peak from Gaussian originals, then calls the function
% DEMG1stDerivativeManualSymm to perform the symmetrization, using a
% three-level plus-and-minus bracketing technique to help you to determine
% the best values of the weighting factors.
% Requires the following functions: expgaussian2, ExpBroaden, fastsmooth,
% halfwidth, peakfit, autopeaks.
format compact
format short g
A=100; % Peak amplitude of peak 1
s=100; % Standard deviation of underlying Gausssian of EMG
lambda1=.01; % first exponential modification
lambda2=.003; % second exponential modification
mu=500; % x-axis location of first peak
separation=900; % Separation between peaks
A2=50; % Relative height of seconf peak
deltat=5; % x-axis increment
Noise=.01; % Random white noise
smoothwidth=5; % Smoothing applied to signal 
t=0:deltat:3500; % Create time axis

% Set the center values of the brackets for L2 and L1
L2=333; % Adjust this first (set L1 to zero, then adjust L2)
L1=100;  % Adjust this after L2 is optimized

Delta=2; % Produces +/-1% brackets

% Construction of exponentially modified Gaussian (EMG)
tau1=1/lambda1;
EMG=A.*expgaussian2(t,mu,s,tau1)+A2.*expgaussian2(t,mu+separation,s,tau1)+Noise.*randn(size(t));
EMG=fastsmooth(EMG,smoothwidth,3);
% Determine true areas under peaks
TrueArea1=trapz(t,A.*expgaussian2(t,mu,s,tau1));
TrueArea2=trapz(t,0.5.*A.*expgaussian2(t,mu+separation,s,tau1));

% Application of second stage of exponential broadening to the EMG
tau2=1/lambda2;
DEMG=ExpBroaden(EMG',-tau2./deltat)'; % DEMG=Doubly Exponential Mod. Gaussian

% Call DEMSymm function
rG=DEMSymm(t,DEMG,L1,L2,Delta);
% Repeat this function with improved values of L1guess and L2guess until
% the middle ones are best (adjust L2 first, then L1)

% disp('        DEMP      Recovered')
% Halfwidth1a=halfwidth(t,DEMG,mu);
% Halfwidth2a=halfwidth(t,rG,mu);
% Halfwidth1b=halfwidth(t,DEMG,mu+separation);
% Halfwidth2b=halfwidth(t,rG,mu+separation);
% HalfWidths=[Halfwidth1a Halfwidth2a]
% HalfWidths=[Halfwidth1b Halfwidth2b]

disp('Peak area measurement of recovered peaks by perpendicular drop method')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
autopeaksresults=autopeaks(t,rG);
disp(autopeaksresults)
disp(' ')
% If the peakfit function is available, see if recovered peak is Gaussian
try
figure(3)
disp('    Fit Gaussians to Recovered peaks')

[FitResults,GOF]=peakfit([t;rG],0,0,2,1,0,10);
disp('          Peak    Position       Height      Halfwidth      Area ')
disp(FitResults)
disp('  % fit error           R2')
disp(GOF)
subplot(2,1,1);title('Recovered peak fit to Gaussians')
catch
end

disp('% area accuracy')
disp(100.*(TrueArea1-autopeaksresults(1,5))./TrueArea1)
disp(100.*(TrueArea2-autopeaksresults(2,5))./TrueArea2)


disp(' ')
