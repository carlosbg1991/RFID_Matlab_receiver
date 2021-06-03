% A script that simulates the calibration curve of three overlapping
% Gaussian peaks. Even-derivative sharpening (the red line in the signal
% plots)is used to improve the resolution of the peaks to allow
% perpendicular drop area measurement. A straight line is fit to the
% calibration curve and the R2 is calculated, in order to demonstrate (1)
% the linearily of the response, and (2) in independence on the overlapping
% peak 1. Must have gaussian.m, derivxy.m, autopeaks.m, val2ind.m,
% halfwidth.m, and plotit.m in the path.
format compact
format short g
clear
% User-controllable variables
scale=1; % x-axis scaling factor
dx=0.001*scale;  % Sampling interval
w=2*scale; % FWHM (full width at half maximum). <<< Yalue of 2 gives Rs=0.55
pos1=16*scale; % Position (retention time) of peak 1 
pos2=18.2*scale; % Position (retention time) of peak 2 
pos3=20.4*scale;% Position (retention time) of peak 3 
maxamp=1; % Minimum amplitude of sample peaks <<<
minamp=0.2; %  Maximum amplitude of sample peaks <<<
factor=32; % Sharpening control for calculating k2 and k4 on lines 33, 34
MaxTrials=4; % <<< Number of standards on the calibration curve. <<<
% Peak detection parameters
SlopeThreshold=.00000000001;
AmpThreshold=.1;
smoothwidth=1;
peakgroup=5;
    
% Calculations
x=11*scale:dx:28*scale;  % time axis
Resolution1=(pos2-pos1)/(2*w);
Resolution2=(pos3-pos2)/(2*w);
k2=1.7*((w).^2)/factor; 
k4=1.7*((w).^4)/(28*factor); 
disp(['K2 = ' num2str((1.7*(w).^2)/factor) ] )
disp(['K4 = ' num2str(1.7*((w).^4)/(28*factor)) ] )
figure(1)   
clf
disp('Peak amplitudes and measured areas of the standards')
disp('          amp1        amp2       amp3         Area 1       Area 2       Area3 ')
for trial=1:MaxTrials
    amp1(trial)=minamp+maxamp.*rand(); % Amplitude range of peak 1
    amp2(trial)=minamp+maxamp.*rand(); % Amplitude range of peak 2
    amp3(trial)=minamp+maxamp.*rand(); % Amplitude range of peak 2
    y1=amp1(trial).*gaussian(x,pos1,w); %  isolated peak 1
    y2=amp2(trial).*gaussian(x,pos2,w); %  isolated peak 2
    y3=amp3(trial).*gaussian(x,pos3,w); %  isolated peak 3
    y=y1+y2+y3; % Measured signal is the sum of the 3 peaks
    % Actual area of each isolated peak
    TrueArea1(trial)=trapz(x,y1);
    TrueArea2(trial)=trapz(x,y2);
    TrueArea3(trial)=trapz(x,y3);
    % Predicted area based on measured value of y at pos1 and pos2, assuming Gaussian.
    % Used only of the perpendicular drop method fails due to insufficient
    % valley
    AproxArea1(trial)=1.0645*(y(val2ind(x,pos1)))*w;
    AproxArea2(trial)=1.0645*(y(val2ind(x,pos2)))*w;
    AproxArea3(trial)=1.0645*(y(val2ind(x,pos3)))*w;
    % Derivative and Sharpening Calculations
    d=derivxy(x,y); % First derivative term
    d2=derivxy(x,d);% Second derivative term
    d3=derivxy(x,d2);% Third derivative term
    d4=derivxy(x,d3); % Fourth derivative term
    ey=y-k2.*d2+k4.*d4; % ey=sharpened signal
    % Prepare subplot array for showing all "MaxTrials" signals
    side=round(sqrt(MaxTrials)+.4); % Number of subplots on a side of square array
    subplot(side,side,trial)
    plot(x,y,x,ey)
    title(['Height ratios = ' num2str([amp1(trial)./amp1(trial) amp2(trial)./amp1(trial) amp3(trial)./amp1(trial)])]);
    % axis([min(x) max(x) 0 1.6])
    
    try
        APresults=autopeaks(x,ey,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,3);
        PerpDrop=APresults(:,5);
        area1(trial)=PerpDrop(1);
        area2(trial)=PerpDrop(2);
        area3(trial)=PerpDrop(3);
    catch
        disp('autopeaks failed, PD area replaced by estimate from 1.0645*h*w for Gaussian')
        sizer=size(APresults);PeaksDetected=sizer(1)
        area1(trial)=AproxArea1(trial);
        area2(trial)=AproxArea2(trial);
        area3(trial)=AproxArea3(trial);
    end
    disp([amp1(trial) amp2(trial) amp3(trial) area1(trial) area2(trial) area3(trial)] )
end
figure(2)
[coef1,RSquared,StdDevs1]=plotit(amp1,area1,1,'o');
title(['Calibration curve.   Rs1 = ' num2str(Resolution1) '    Rs2 = ' num2str(Resolution2) '    Factor = ' num2str(factor) ])
ylabel('Peak area')
xlabel('amplitude of peak 1')
% Rsquared1=RSquared
% Use the standards as unknowns and analyze them using this calibration curve.
PredAmp1=(area1-coef1(2))./coef1(1); % Predicted amp1 from calibration curve
PercentErrors1=100.*(amp1-PredAmp1)./amp1;
AveragePercentError1=mean(abs(PercentErrors1));

figure(3)
[coef2,RSquared,StdDevs2]=plotit(amp2,area2,1,'o');
title(['Calibration curve.   Rs1 = ' num2str(Resolution1) '    Rs2 = ' num2str(Resolution2) '    Factor = ' num2str(factor) ])
ylabel('Peak area')
xlabel('amplitude of peak 2')
% Rsquared2=RSquared
% Use the standards as unknowns and analyze them using this calibration curve.
PredAmp2=(area2-coef2(2))./coef2(1); % Predicted amp2 from calibration curve
PercentErrors2=100.*(amp2-PredAmp2)./amp2;
AveragePercentError2=mean(abs(PercentErrors2));

figure(4)
[coef3,RSquared,StdDevs3]=plotit(amp3,area3,1,'o');
title(['Calibration curve.   Rs1 = ' num2str(Resolution1) '    Rs = ' num2str(Resolution2) '    Factor = ' num2str(factor) ])
ylabel('Peak area')
xlabel('amplitude of peak 3')
% Rsquared3=RSquared
% Use the standards as unknowns and analyze them using this calibration curve.
PredAmp3=(area3-coef3(2))./coef3(1); % Predicted amp3 from calibration curve
PercentErrors3=100.*(amp3-PredAmp3)./amp3;
AveragePercentError3=mean(abs(PercentErrors3));

disp('Average % errors for analyzing the standards as unknowns')
disp('    Component 1  Component 2  Component 3')
disp([AveragePercentError1 AveragePercentError2 AveragePercentError3])

disp('Maximum % errors for analyzing the standards as unknowns')
disp('    Component 1  Component 2  Component 3')
disp([max(abs(PercentErrors1)) max(abs(PercentErrors2)) max(abs(PercentErrors3))])
disp(' ')