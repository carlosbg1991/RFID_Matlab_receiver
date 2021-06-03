% Demonstration of peak measurement with very low signal-to-noise ratios.
% Change MaxPeakHeight (line 8) to make the problem harder or easier.
% Compares smoothing, findpeaksG, peakfit, and classical lest squares (cls2.m).
increment=0.1;% Change increment to change the number of points in signal.
x=-100:increment:300; 
Noise=1; % Leave this at 1 so the SNR equals the peak height on the plots.
Baseline=1; % Change if you wish (AmpThreshold will adjust automatically)
MaxPeakHeight=2; % Change if you wish. SNR range will be from zero to this value.
PeakPosition=100; % Change to test effect of shifts in peak position
PeakWidth=100; % Change to test effect of shifts in peak width
for trial=1:41,
    PeakHeight(trial)=MaxPeakHeight./40*(trial-1);
    SignalToNoiseRatio=PeakHeight./Noise;
    % For a Lorentzian peak, change 'gaussian' to 'lorentzian' in line 16,
    % 'findpeaksG' to 'findpeaksL' in line 28 and PeakShape (line 22) to 2.
    y=Baseline+PeakHeight(trial).*gaussian(x,PeakPosition,PeakWidth)+Noise.*randn(size(x));
    SlopeThreshold=0.2*(PeakWidth./increment)^-2; 
    AmpThreshold=Baseline+MaxPeakHeight./50;
    SmoothWidth=PeakWidth./increment; 
    FitWidth=PeakWidth./increment; 
    SmoothMode=3;
    PeakShape=1; % Test the effect of a shape mismatch by changing this.
    NumPeaks=1;
    BaselineMode=3; % BaselineMode=3 automatically corrects for flat baseline
    smoothedy=fastsmooth(y,SmoothWidth,3,1);
    SmoothPeakMax(trial)=max(smoothedy)-min(smoothedy);
    SmoothPeakMaxPosition(trial)=x(val2ind(smoothedy,max(smoothedy)));
    FindPeaks=findpeaksG(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,SmoothMode);
    % FindPeaks=findpeaksfit(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,SmoothMode,PeakShape,0,1,1,0,1);
    if FindPeaks==0,
%         FindPeaksHeight(trial)=0;
%         FindPeaksPosition(trial)=0;
    else
        FindPeaksHeight(trial)=FindPeaks(val2ind(FindPeaks(:,2),PeakPosition),3);
        FindPeaksPosition(trial)=FindPeaks(val2ind(FindPeaks(:,2),PeakPosition),2);
    end
    PeakFit=peakfit([x;y],0,0,NumPeaks,PeakShape,0,1,0,BaselineMode);
    FitHeight(trial)=PeakFit(3);
    FitPosition(trial)=PeakFit(2);
    PowerSpectrum=PlotFrequencySpectrum(x',y',0,0,1);
    CLSresults=cls2(x,y,1,1,100,100,0);
    CLSheight(trial)=CLSresults(2);
    PS(trial)=PowerSpectrum(2,2); % First harmonic
end
SmoothPeakMaxError=100.*mean((abs(SmoothPeakMax-PeakHeight))/MaxPeakHeight);
FindPeaksHeightError=100.*mean((abs(FindPeaksHeight-PeakHeight))/MaxPeakHeight);
FitHeightError=100.*mean((abs(FitHeight-PeakHeight))/MaxPeakHeight);
SmoothPeakMaxPositionError=100.*std(SmoothPeakMaxPosition)./mean(SmoothPeakMaxPosition);
FindPeaksPositionError=100.*std(FindPeaksPosition)./mean(FindPeaksPosition);
FitPositionError= 100.*std(FitPosition)./mean(FitPosition);
CLSHeightError=100.*mean((abs(CLSheight-PeakHeight))/MaxPeakHeight);
disp('-------------------------------------------------------')
disp(['Number of points in half-width of peak: '  num2str(PeakWidth./increment)  ])
disp('Method         Height Error     Position Error')
disp(['Smoothed peak    ' num2str(SmoothPeakMaxError) '%        ' num2str(SmoothPeakMaxPositionError) '%' ] )
disp(['findpeaksG.m     ' num2str(FindPeaksHeightError) '%        ' num2str(FindPeaksPositionError) '%' ] )
disp(['peakfit.m        ' num2str(FitHeightError) '%         '  num2str(FitPositionError) '%' ] )
disp(['cls2.m           ' num2str(CLSHeightError) '%         ' ] )

clf
figure(1) % Peak height measurements compared in Figure 1.
subplot(2,2,2)
plotit(PeakHeight,FindPeaksHeight,1);
xlabel('Actual peak height (also signal-to-noise ratio)')
ylabel('Measured peak height')
title('findpeaksG function')
subplot(2,2,1)
plotit(PeakHeight,SmoothPeakMax,1);
xlabel('Actual peak height (also signal-to-noise ratio)')
ylabel('Measured peak height')
title('Maximum of smoothed peak')
subplot(2,2,3)
plotit(PeakHeight,FitHeight,1);
xlabel('Actual peak height (also signal-to-noise ratio)')
ylabel('Measured peak height')
title('Unconstrained Gaussian least-squares peak fit')
subplot(2,2,4)
plotit(PeakHeight,CLSheight,1);
xlabel('Actual peak height (also signal-to-noise ratio)')
ylabel('Measured peak height')
title('Gaussian least squares fit: position=100, width=100)')
figure(2)
plotit(PeakHeight,sqrt(PS),1);
xlabel('Actual peak height ')
ylabel('Measured fundamental frequency amplitude')
title('Fundamental frequency amplitude')
AXIS([0 2 0 1700])
figure(1)