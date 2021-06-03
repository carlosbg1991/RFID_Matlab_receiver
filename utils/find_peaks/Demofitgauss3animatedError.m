% Demos iterative gaussian fit to three computer-generated noisy gaussian
% peaks. Variable parameters are height, position, and width of each peak.
% Requires the gaussian.m, fitgauss3animatedError.m, and fminsearch.m
% functions to be in the Matlab path. Thise version adds a plot of fitting
% error vs iteration nunber. T. C. O'Haver, June 2019

format compact
global c NumTrials TrialError
warning off
% Create a three gaussian peaks and add random noise
k=10; % sets interval between x-values; smaller k = more points in signal
NumTrials=0;
TrialError=0;
t=-0:k:600;
noise=0;
Height1=100;
Position1=200; 
Width1=90;  
Height2=50;
Position2=300; 
Width2=100;  
Height3=30;
Position3=400; 
Width3=110;  
peak1=[Position1 Height1 Width1]
peak2=[Position2 Height2 Width2]
peak3=[Position3 Height3 Width3]
y = Height1.*gaussian(t,Position1,Width1) + Height2.*gaussian(t,Position2,Width2)+ Height3.*gaussian(t,Position3,Width3) + noise.*randn(size(t));

% Perform an iterative fit using the FMINSEARCH function
% Estimate start values from the independent axis (t values) alone.
mint=min(t);
maxt=max(t);
dt=maxt-mint;
p1=150;
p2=300;
p3=500;
w1=dt/5;
w2=dt/5;
w3=dt/5;
% Starting point "start" has the form [position1 width1 position2 width2]
start=[p1 w1 p2 w2 p3 w3]  % a reasonable guess
% start=[p1 w1 p1 w1 p1 w1]  % really bad guess
% start=[min(t) max(t)./10 max(t)./2 max(t)./10 max(t) max(t)./10]  %  really bad guess
options = optimset('TolX',0.1,'MaxFunEval',600);  % Determines how close the model must fit the data
% "tic" and "toc" measure how long the iterative fit takes.
% tic
parameter=fminsearch(@(lambda)(fitgauss3animatedError(lambda,t,y)),start);
% toc
disp('      Height      Position      Width')
MeasuredPeakMatrix=[parameter(1) c(1) parameter(2);parameter(3) c(2) parameter(4); parameter(5) c(3) parameter(6)] % Peak height 2 is returned in the global variable c(2).

% Compute a model and plot it (blue line) along with 
% the original data (red points)
Peak1=c(1).*gaussian(t,parameter(1),parameter(2));
Peak2=c(2).*gaussian(t,parameter(3),parameter(4));
Peak3=c(3).*gaussian(t,parameter(5),parameter(6));
model=Peak1+Peak2+Peak3;
clf
subplot(1,2,1)
plot(t,y,'ro',t,Peak1,'c',t,Peak2,'m',t,Peak3,'g',t,model,'b')

subplot(1,2,2)
semilogy(1:NumTrials,TrialError)
grid
xlabel('Iteration Number')
ylabel('Error')
NumTrials
title('Progress toward minimizing the fitting error')
disp('Percent parameter errors:')
disp('      Height      Position      Width')
disp(100.*(peak1-MeasuredPeakMatrix(1,:))./peak1);
disp(100.*(peak2-MeasuredPeakMatrix(2,:))./peak2);
disp(100.*(peak3-MeasuredPeakMatrix(3,:))./peak3);
disp(' ')
