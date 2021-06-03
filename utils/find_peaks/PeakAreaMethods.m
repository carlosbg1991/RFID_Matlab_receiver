% Area measurements of 5 overlapping peaks with unequal heights but equal
% areas. Areas are measured by Gaussian estimation (GE), perpendicular drop
% (PD), and iterative curve fitting (Fa), which are taken as the correct
% values. (All three of these methods are available in the interactive peak
% finder iPeak.m). Try different peak widths (line 19) and shapes (line 18).
%
% Conclusions: 
% 1. Area measurements are more accurate when peaks are more
% separate and less overlapped. 
% 2. For slightly overlapped Gaussian peaks (PeakShape=1) both the Gaussian
% estimation and perpendicular drop methods work well, but for Logistic and
% Lorentzian shapes, or for heavily overlapped Gaussians, perpendicular
% drop is much better.
%
% Requires the modelpeaks, autofindpeaks, autopeaks, ipeak, and peakfit
% functions to be in the Matlab path.
x=1:.1:70;
PeakShape=1; % <<< Gaussian=1, Logistic distribution=3, Lorentzian=2, etc
width=7; % <<< larger values = more overlap and harder measurement.
% calculations
shapes=[1 1 1 1 1].*PeakShape;
heights=[1 1.1 1.2 1.3 1.4];
positions=[55 45 35 25 15];
widths=width./heights; % Widths inversely proportional to heights to keep areas constant
extras=[0 0 0 0 0];
start=reshape([positions' widths']',1,10);
y=modelpeaks2(x, shapes, heights, positions, widths, extras);

ipeak(x,y,length(shapes));

P=autofindpeaks(x,y,7);
Pa=P(:,5);

M=autopeaks(x,y,7);
Ma=M(:,5);

[FitResults,FitError]=peakfit([x;y],0,0,5,PeakShape,0,1,start,0,0,0);
Fa=FitResults(:,5);
 
PercentErrorsGE=100*(Pa-Fa)./Fa; 
PercentErrorsPD=100*(Ma-Fa)./Fa;

GaussianEstimationErrors=[[1:5]' PercentErrorsGE]
PerpendicularDropErrors=[[1:5]' PercentErrorsPD]