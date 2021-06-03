% Script to test segmented peak finder findpeaksSG. Requires findpeaksSG.m
% and gaussian.m (or other peak shape) function. Creates 100 noisy signals
% with 3 peaks with different widths, then measures the peak positions
% heights and widths, and computes the percent relative standard deviation
% of 100 measurements. The segmented peak finder findpeaksSG, with
% 3-segment peak detection parameters, reliably detects amd measures all
% three peaks. In contrast, findpeaksG tuned to the middle peak (using line
% 26 instead of line 25) measures the first and last peaks poorly. You can
% change any of the variables in lines 10-18.
clear
clf
p=[200 400 800]; % Peak positions of the 3 peaks.
h=[1 1 1]; % Peak heights of the 3 peaks.
w=[10 50 200]; % Peak widths of the 3 peaks.
Noise=.04; % Standard deviation of the random noise
% Define the peak detection parameters for each segment
SlopeThresholds=[.001 .00001 .000001];
AmpThresholds=[.1 .1 .1];
smoothwidths=[10 50 80];
peakgroups=[8 50 250];
NumTrials=100; % Number of independent signals measured
x=1:1000;
for trial=1:NumTrials,
    % Create the noisy signal
    y=h(1).*gaussian(x,p(1),w(1))+h(2).*gaussian(x,p(2),w(2))+h(3).*gaussian(x,p(3),w(3))+Noise.*randn(size(x));
    % Detect and measure peaks in the signal
    P=findpeaksSG(x,y,SlopeThresholds,AmpThresholds,smoothwidths,peakgroups,3);
    % P=findpeaksG(x,y,SlopeThresholds(2),AmpThresholds(2),smoothwidths(2),peakgroups(2),3);
    sizeP=size(P);
    NumPeaks(trial)=sizeP(1);
    if trial==1; % Plot the first signal
        plot(x,y,'c'),
    end
    hold on
    for peak=1:NumPeaks(trial),
        text(P(peak,2),P(peak,3),num2str(peak)) % Number the peaks
        % Record the measured peak parameters of each peak in each trial
        PeakPositions(trial,peak)=P(peak,2);
        PeakHeights(trial,peak)=P(peak,3);
        PeakWidths(trial,peak)=P(peak,4);
        PeakAreas(trial,peak)=P(peak,5);
    end
end
hold off
xlabel('Each detected peak position and height is marked with peak number 1, 2 or 3')
ylabel('Y')
title('TestPrecisionFindpeaskSG.m: 100 repeat measurements of noisy signal with 3 peaks')
disp(' ')
disp(['Average number of peaks: ' num2str(mean(NumPeaks)) ])
disp('Percent Relative Standard Deviations (RSD):')
disp('     Peak 1      Peak 2      Peak 3       Peak 4')
RSDPeakPositions=100.*std(PeakPositions)./mean(PeakPositions)
RSDPeakHeights=100.*std(PeakHeights)./mean(PeakHeights)
RSDPeakWidths=100.*std(PeakWidths)./mean(PeakWidths)
RSDPeakAreas=100.*std(PeakAreas)./mean(PeakAreas)