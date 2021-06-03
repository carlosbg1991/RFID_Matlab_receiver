% Demonstration of symmetrization of a double exponentially modified
% Gaussian (EMG) by two successive applicaitons of weighted addition of the
% first derivative, with weighting factors equal to the two taus. Shows
% that the peak width is reduced, the peak made more symmetrical, and the
% area is preserved. This demo scripts creates a double exponential peak
% from a Gaussian original, then calls DEMG1stDerivativeManualSymm to
% perform the symmetrization, using a three-level plus-and-minus bracketing
% technique to help you to determine the best weighting factors.
% Requires the following functions: expgaussian2, ExpBroaden, fastsmooth,
% halfwidth, peakfit, autopeaks.
format compact
format short g
A=100; % Peak amplitude
s=100; % Standard deviation of underlying Gausssian of EMG
lambda1=.01; % first exponential modification
lambda2=.003; % second exponential modification
mu=500; % x-axis location of first peak
deltat=5; % x-axis increment
Noise=.01; % Random white noise
smoothwidth=5; % Smoothing applied to signal 
t=0:deltat:3500; % Create time axis

% Set the center values of the brackets for L2 and L1
L2=333; % Adjust this first (set L1guess to zero)
L1=100;  % Adjust this after L2 is optimized

Delta=1; % Produces +/-1% brackets 

% Construction of exponentially modified Gaussian (EMG)
tau1=1/lambda1;
EMG=A.*expgaussian2(t,mu,s,tau1)+Noise.*randn(size(t));
EMG=fastsmooth(EMG,smoothwidth,3);
% Determine true area under peak
TrueArea=trapz(t,A.*expgaussian2(t,mu,s,tau1));

% Application of second stage of exponential broadening to the EMG
tau2=1/lambda2;
DEMG=ExpBroaden(EMG',-tau2./deltat)'; % DEMG=Doubly Exponential Mod. Gaussian

% Call DEMSymm function
rG=DEMSymm(t,DEMG,L1,L2,Delta);
% Repeat this function with improved values of L1guess and L2guess until
% the middle ones are best (adjust L2 first, then L1)

disp('        DEMP     Recovered')
Halfwidth1=halfwidth(t,DEMG);
Halfwidth3=halfwidth(t,rG);
HalfWidths=[Halfwidth1 Halfwidth3]

Area1=trapz(t,DEMG);
Area3=trapz(t,rG);
Areas=[Area1 Area3];
disp(' ')

% If the peakfit function is available, see if recovered peak is Gaussian
try
figure(3)
disp('    Fit Gaussian to Recovered peak')

[FitResults,GOF]=peakfit([t;rG]);
disp('          Peak    Position       Height      Halfwidth      Area ')
disp(FitResults)
disp('  % fit error           R2')
disp(GOF)
subplot(2,1,1);title('Recovered peak fit to Gaussian')
catch
end
disp('% area accuracy')
disp(100.*(Area3-Area1)./Area1)
disp(' ')

