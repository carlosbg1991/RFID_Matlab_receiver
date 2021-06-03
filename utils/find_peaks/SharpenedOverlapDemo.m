% A script that automatically determines the optimum degree of even-
% derivative sharpening for two overlapping Gaussians that minimizes the
% errors of measuring peak areas by perpedicular drop. Must have
% gaussian.m, derivxy.m, autopeaks.m, val2ind.m, and halfwidth.m in the
% path. Download these from https://terpconnect.umd.edu/~toh/spectrum/
format compact
format short g
clear
scale=1; % x-axis scaling factor
dx=0.001*scale;  % Sampling interval
x=11*scale:dx:25*scale;  % time axis
w=1.5*scale; % FWHM (full width at half maximum) <<<<< CHANGE RESOLUTION HERE
amp1=1; % Amplitude of peak 1
amp2=.5;% Amplitude of peak 2  <<<<< CHANGE AMPLITUDE HERE
pos1=16*scale; % Position (retention time) of peak 1
pos2=18.2*scale; % Position (retention time) of peak 2  
y1=amp1.*gaussian(x,pos1,w);
y2=amp2.*gaussian(x,pos2,w);
y=y1+y2;
% y=y+0.00001*randn(size(y));
% y=fastsmooth(y,1001,3,0);
TrueArea1=trapz(x,y1);
TrueArea2=trapz(x,y2);
Resolution=(pos2-pos1)/(2*w);
trial=0;
for factor=70:-.5:20
    trial=trial+1;
    factorvector(trial)=factor;
    % Estimated values of weighting factors k2 and k4 based on peak width
    k2(trial)=1.7*((w).^2)/factor; % estimate k2 
    k4(trial)=1.7*((w).^4)/(28*factor); % estimate k4
    % Calculations
    d=derivxy(x,y); % First derivative term
    d2=derivxy(x,d);% Second derivative term
    d3=derivxy(x,d2);% Third derivative term
    d4=derivxy(x,d3); % Fourth derivative term
    ey=y-k2(trial).*d2+k4(trial).*d4; % ey=sharpened signal
    % ey=(y-k2(trial).*d2+k4(trial).*d4).*(y/max(y)); % ey=sharpened signal
    TwoWarea1b=sum(y(val2ind(x,pos1-w):(val2ind(x,pos1+w)))).*dx;
    TwoWarea2b=sum(y(val2ind(x,pos2-w):(val2ind(x,pos2+w)))).*dx;
    TwoWarea1=sum(ey(val2ind(x,pos1-w):(val2ind(x,pos1+w)))).*dx;
    TwoWarea2=sum(ey(val2ind(x,pos2-w):(val2ind(x,pos2+w)))).*dx;
    
% Remove commenting below to plot animation of signals with different sharpening values
%     figure(1)
%     clf
%     plot(x,y,x,ey)
%     % axis([10.2 25.8 -.5 2])
%     xlabel(['Blue=Original Gaussians, Resolution, Rs = ' num2str(Resolution) '     factor = ' num2str(factor)  ])
%     ylabel('signal amplitude')
%     title('Effect of sharpening on perpendicular drop area measurements of two overlapping Gaussians')
    
    FWHMy1=halfwidth(x,y,pos1);
    FWHMey1=halfwidth(x,ey,pos1);
    FWHMy2=halfwidth(x,y,pos2);
    FWHMey2=halfwidth(x,ey,pos2);
    yArea=trapz(x,y);
    eyArea=trapz(x,ey);
    PercentChangeInPeak1Width=100.*(FWHMy1-FWHMey1)/FWHMy1;
    PercentAreaError1b=100*(TrueArea1-TwoWarea1b)/TrueArea1;
    PercentAreaError2b=100*(TrueArea2-TwoWarea2b)/TrueArea2;
    PercentAreaError1=100*(TrueArea1-TwoWarea1)/TrueArea1;
    PercentAreaError2=100*(TrueArea2-TwoWarea2)/TrueArea2;
    
    SlopeThreshold=.0000000001;
    AmpThreshold=.01;
    smoothwidth=101;
    peakgroup=51;
    try
        APresultsb=autopeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,3);PerpDropb=APresultsb(:,5);
        PDError1b(trial)=100*(PerpDropb(1)-TrueArea1)/TrueArea1;
        PDError2b(trial)=100*(PerpDropb(2)-TrueArea2)/TrueArea2;
    catch
        disp('perpendicular drop failed')
        PDError1b(trial)=0;
        PDError2b(trial)=0;
    end
    try
        APresults=autopeaks(x,ey,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,3);PerpDrop=APresults(:,5);
        PDError1(trial)=100*(PerpDrop(1)-TrueArea1)/TrueArea1;
        PDError2(trial)=100*(PerpDrop(2)-TrueArea2)/TrueArea2;
    catch
        disp('perpendicular drop  failed')
        PDError1(trial)=0;
        PDError2(trial)=0;
    end
    valley(trial)=min(ey(val2ind(x,pos1):val2ind(x,pos2)));
end
% Find factors at which minimum errors occurs
factor1=factorvector(val2ind(PDError1,min(abs(PDError1))));
factor2=factorvector(val2ind(PDError2,min(abs(PDError2))));
AverageError=(abs(PDError1)+abs(PDError2))/2;
AverageFactor=factorvector(val2ind(AverageError,min(AverageError)));
factor=AverageFactor;
% Perform calculation at this value
k2=1.7*((w).^2)/factor;
k4=1.7*((w).^4)/(28*factor); 
ey=y-k2.*d2+k4.*d4; % ey=sharpened signal 

% Plot signals with optimum sharpening for minimum errors
figure(1)
clf
plot(x,y,x,ey)
xlabel(['Blue=Original Gaussians, Resolution, Rs = ' num2str(Resolution) '     factor = ' num2str(factor) '    Red: Sharpened peaks' ])
ylabel('signal amplitude')
title('Optimum sharpening for minimum errors of perpendicular drop area measurements')
trial=val2ind(PDError2,min(abs(PDError2)));
text(pos1,amp1,['    Area error = ' num2str(PDError1b(trial))  '%'] )
text(pos1,max(ey),['    Area error = ' num2str(PDError1(trial))  '%'] )
text(pos2,amp2,['    Area error = ' num2str(PDError2b(trial)) '%'] )
text(pos2,ey(val2ind(x,pos2)),['    Area error = ' num2str(PDError2(trial)) '%'] )

figure(2)
plot(factorvector,abs(PDError1),factorvector,abs(PDError2),factorvector,abs(valley)/max(abs(valley)),factorvector,AverageError,'k')
title(['Blue: peak 1.  Red: peak 2.  Black: average   Yellow: valley.   Rs = ' num2str(Resolution) '     Amp. ratio = ' num2str(amp1/amp2) ])
ylabel('Absolute percent area error')
xlabel('Weighting factor')
disp(['Optimum factor = ' num2str(factor) ] )
disp(['Optimum K2 = ' num2str((1.7*(w).^2)/factor) ] )
disp(['Optimum K4 = ' num2str(1.7*((w).^4)/(28*factor)) ] )