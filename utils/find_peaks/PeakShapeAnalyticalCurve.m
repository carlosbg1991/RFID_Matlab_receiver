% Demonstration of a quantitative analysis application of curve fitting,
% where the peak height measured by curve fitting is used only to determine
% the concentration of the substance that created the peak by constructing
% a calibration curve. This demo shows that the peak shape is relatively
% unimportant, for a single isolated peak whose shape is constant
% and independent of concentration. If the wrong model shape is used, the
% peak heights measured by curve fitting will be inaccurate, but the error
% will be exactly the same for the unknown samples and the known
% calibration standards, so the error will cancel out and the measured
% concentrations will be accurate. Still, it's useful to use as accurate a
% model peak shape as possible, because the goodness of fit (percent
% fitting error or the R2) can be used as an indicator if something
% unexpected goes wrong during the analysis (such as an increase in the
% noise or the appearance of an interfering peak from a foreign substance).
t=[1:5:1000]; % Defines the "number of wavelengths"
Noise=.1;
FitPeakShape=1; % Try the "wrong" peakshape here (try 2,3,21,41)
sensitivity=10; % proportionality constant between concentration and peak height
% The number of measurements is defined by the length of amp
conc=[.1 .1 1 1 2 2 3 3 4 4 5 5];  % concentrations of standards
pos=500;   % Position of peak  on the x axis
wid=200;   % Width (full width at half maximum) of the peak
NumPeaks=1;NumTrials=1;extra=0;start=0;autozero=1;fixedparameters=0;plots=1;bipolar=0;minwidth=.1;DELTA=1;clipheight=0;
% A = matrix containing one unit-amplitude peak in each of its rows
A = zeros(length(pos),length(t));
amp=(conc*sensitivity); % amp is directly proportional to conc.
% amp=(conc*sensitivity).^.9; % Nonlinear relationship between conc and amp
figure(1)
for k=1:length(amp)
 A(k,:)=t/100+amp(k).*gaussian(t,pos,wid)+Noise.*randn(size(t)); % You can use any peak function here
[results(k,:),R2(k,:)]=peakfit([t;A(k,:)],0,0,NumPeaks,FitPeakShape,extra,NumTrials,start,autozero);
% [results(k,:),R2(k,:)]=peakfit([t;A(k,:)],0,0,NumPeaks,FitPeakShape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight);
end
subplot(2,1,1);title('Curve fitting result for the last standard (using peakfit.m)')
MeanR2=mean(R2(:,2));
figure(2)
clf
plot(t,A)
xlabel('x')
ylabel('Signal Amplitude')
title('Superimposed signal peaks of all the standards')
figure(3)
clg
coef=plotit(conc,results(:,3),1);
title('Calibration curve based on peak heights measured by curve fitting')
MeasuredConc=(results(:,3) - coef(2))./coef(1);
xlabel('Concentrations of Standards')
ylabel('Measured peak heights')
figure(4)
plotit(conc,MeasuredConc,1);
title('Correlation of actual and measured concentrations, based on peak height')
xlabel('Actual concentrations')
ylabel('Measured concentrations')
figure(5)
clg
coef=plotit(conc,results(:,5),1);
title('Calibration curve based on peak areas measured by curve fitting')
ylabel('Measured peak areas')
xlabel('Concentrations of Standards')
% Using the peak area readings of the standards as unknown readings
MeasuredConc=(results(:,5) - coef(2))./coef(1);
figure(6)
plotit(conc,MeasuredConc,1);
title('Correlation of actual and measured concentrations, based on peak area')
xlabel('Actual concentrations')
ylabel('Measured concentrations')