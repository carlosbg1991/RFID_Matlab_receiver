% Demos iterative gaussian fit to two computer-generated noisy gaussian
% peaks. Variable parameters are height, position, and width of both peaks.
% Requires the gaussian, fitgauss2animated, and fminsearch functions to be
% in the Matlab path. T. C. O'Haver, May 2019
format compact
global c NumTrials TrialError
warning off
% Create a two gaussian peaks and add random noise
k=10; % sets interval between x-values; smaller k = more points in signal
NumTrials=0;
TrialError=0;
t=-0:k:500;
noise=1;
Height1=100;
Position1=200; 
Width1=200;  
Height2=50;
Position2=400; 
Width2=100;  
Peak1=[Height1 Position1 Width1]
Peak2=[Height2 Position2 Width2]
y = Height1.*gaussian(t,Position1,Width1) + Height2.*gaussian(t,Position2,Width2) + noise.*randn(size(t));

% Perform an iterative fit using the FMINSEARCH function
% Estimate start values from the independent axis (t values) alone.
mint=min(t);
maxt=max(t);
dt=maxt-mint;
p1=mint+dt/3;
p2=mint+dt/2;
w1=dt/5;
w2=dt/5;
% Starting point "start" has the form [position1 width1 position2 width2]
start=[p1 w1 p2 w2]  % automatically calculcated start (not very good)
% start=[150 150 350 150] % Pretty good guess
options = optimset('TolX',0.1);  % Determines how close the model must fit the data
% "tic" and "toc" measure how long the iterative fit takes.
% tic
parameter=fminsearch(@(lambda)(fitgauss2animated(lambda,t,y)),start);
% toc
MeasuredPeak1=[c(1) parameter(1) parameter(2)] % Peak height 1 is returned in the global variable c(1).
MeasuredPeak2=[c(2) parameter(3) parameter(4)] % Peak height 2 is returned in the global variable c(2).

% Compute a model and plot it (blue line) along with 
% the original data (red points)
peak1=c(1).*gaussian(t,parameter(1),parameter(2));
peak2=c(2).*gaussian(t,parameter(3),parameter(4));
model=peak1+peak2;
figure(1)
plot(t,y,'ro',t,peak1,'c',t,peak2,'m',t,model,'b')
disp(' ')
disp('    Height   Position   Width')
Peak1PercentErrors=100.*(Peak1-MeasuredPeak1)./Peak1
Peak2PercentErrors=100.*(Peak2-MeasuredPeak2)./Peak2
disp(' ')