% Comparison of peak height and peak area measurements for a single
% isolated peak with variable shape (caused by variable exponential
% broadening), position, width, and background. You can change any of the
% constants in lines 7-16. Requires the modelpeaks.m, plotit.m, and peakfit.m functions.
clf
clear
deltaX=.1; %  x-axis interval between adjacent points
x=5:deltaX:50; % Creates time axis
standards=[.1 .1 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9]; % True peak heights
% standards=5.*ones(1,20); % Constant peak heights
pos=25; % Position of peak
wid=10; % Width of peak
noise=.1; % Random noise added to signal (0 to 1)
baseline=0; % Average baseline
vba=0; % Varibility of baseline from peak to peak (0 to 2)
vbr=0; % Varibility of peak shape from peak to peak (0 to 2)
for run=1:length(standards)
    tc=3+vbr*35*rand(); % tc = time constant for exponential broadening
    % Create simulated peak - Gaussian with variable exponential broadening
    % plus variable baseline and random noise
    y=baseline*vba*rand()+modelpeaks(x,1,5,standards(run),pos,wid,tc)+noise.*randn(size(x));
    y=real(y.^1); % Simulates detector non-linearity when exponent is not 1.0
    y=y+baseline.*cumsum(y)./(20*wid./deltaX); % Adds a small gradually increasing slope to baseline
    A(run,:)=y; % Computes matrix of peaks for plotting on Figure(1)
    background=(mean(y(1:20)))+mean(y(length(y)-20:length(y)))./2;
    PeakMaximum(run)=max((y)-background); % Subtraction of min(y) helps to correct for baseline
    TotalArea(run)=sum(y-background).*deltaX;
    figure(6)
    [FitResults,Error]=peakfit([x;y],0,0,1,1,0,1,0,1,0,0); % Gaussian fit
    % [FitResults,Error]=peakfit([x;y],0,0,1,31,0,5,0,1,0,0);% ExpGaussian fit (variable)
    PeakHeight(run)=FitResults(3);
    PeakWidth(run)=FitResults(4);
    PeakArea(run)=FitResults(5);
    FitError(run)=Error(1);
%   [w,r]=halfwidth(x,y); % Uncomment these 3 lines to follow peak symmetry
%   FWHM(run)=w;;
%   ratio(run)=r;
    plot(x,y)
    axis([min(x) max(x) 0 max(standards)])
end
clf
figure(1)
plot(x,A) % Shows all peaks plotted together
xlabel('time')
ylabel('Signal')
title('Simulated peaks')
figure(2)
plotit(standards,PeakMaximum,1);
axis([min(standards) max(standards) 0 max(PeakMaximum)])
xlabel('standard')
ylabel('Peak Maximum')
title('Simple Peak Maximum')
figure(3)
plotit(standards,TotalArea,1);
axis([min(standards) max(standards) 0 max(TotalArea)])
xlabel('standard')
ylabel('Total Area')
title('Peak Area by Integration (sum(y)*dx)')
 figure(4)
plotit(standards,PeakHeight,1);
axis([min(standards) max(standards) 0 max(PeakHeight)])
xlabel('standard')
ylabel('Peak Height')
title('Peak Height by Curve Fitting')
figure(5)
plotit(standards,PeakArea,1);
axis([min(standards) max(standards) 0 max(PeakArea)])
xlabel('standard')
ylabel('Peak Area')
title('Peak Area by Curve Fitting')

PercentFittingError=mean(FitError);