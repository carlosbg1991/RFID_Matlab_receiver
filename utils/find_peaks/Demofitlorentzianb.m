% Demos iterative lorentzian fit to two computer-generated noisy lorentzian peaks.
% Variable parameters are height, position, width of both peaks. Baseline
% and peak amplitudes returned in global variable c.
% Uses lorentzian and fitlorentzianb functions.
% T. C. O'Haver, May 2006
format compact
global c

% Create a two lorentzian peaks and add random noise
k=2; % sets interval between x-values; smaller k = more points in signal
t=[-400:k:700];
noise=5;
Height1=100;
Position1=100; 
Width1=100;  
Height2=100;
Position2=200; 
Width2=200;  
Baseline=1;
y = Baseline+Height1.*lorentzian(t,Position1,Width1) + Height2.*lorentzian(t,Position2,Width2) + noise.*randn(size(t));

% Perform an iterative fit using the FMINSEARCH function
start=[length(t)/3 length(t)/2 length(t)/1.5 length(t)/2];    % start now has guesses for BOTH peaks
options = optimset('TolX',0.1);  % Determines how close the model must fit the data
% "tic" and "toc" measure how long the iterative fit takes.
tic
parameter=fminsearch(@(lambda)(fitlorentzianb(lambda,t,y)),start);

toc
Baseline=c(1)
Measured1=[c(2) parameter(1) parameter(2)] % Peak height 1 is returned in the global variable c(1).
Measured2=[c(3) parameter(3) parameter(4)] % Peak height 2 is returned in the global variable c(2).

% Compute a model and plot it (blue line) along with 
% the original data (red points)
peak1=c(2).*lorentzian(t,parameter(1),parameter(2));
peak2=c(3).*lorentzian(t,parameter(3),parameter(4));
model=Baseline+peak1+peak2;
plot(t,y,'r.',t,peak1,'c',t,peak2,'m',t,model,'b',t,Baseline*ones(size(t)))