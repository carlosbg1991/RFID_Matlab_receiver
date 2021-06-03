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
x=10*scale:dx:25*scale;  % time axis
Rs=1; % Chromatographic Resolution <<<<< CHANGE RESOLUTION HERE
factor=40; % <<<<< CHANGE SHARPENING FACTOR HERE
amp1=1; % Amplitude of peak 1
% amp2= is changed in loop, line 24
pos1=16*scale; % Position (retention time) of peak 1
pos2=19*scale; % Position (retention time) of peak 2  
FWHM=((pos2-pos1).*(2*sqrt(2*log(2)))./(4*Rs))*scale; % FWHM (full width at half maximum)
sigma=FWHM/(2*sqrt(2*log(2)));
y1=amp1.*gaussian(x,pos1,FWHM);
TrueArea1=trapz(x,y1);
trial=0;
k2=1.7*((FWHM).^2)/factor;
k4=1.7*((FWHM).^4)/(28*factor); 
for amp2=2:-.02:.05
    trial=trial+1;
    y2=amp2.*gaussian(x,pos2,FWHM);
    TrueArea2=trapz(x,y2);
    y=y1+y2;
    % y=y+0.00001*randn(size(y));
    % y=fastsmooth(y,1001,3,0);
    amp2vector(trial)=amp2;

    % Calculations
    d=derivxy(x,y); % First derivative term
    d2=derivxy(x,d);% Second derivative term
    d3=derivxy(x,d2);% Third derivative term
    d4=derivxy(x,d3); % Fourth derivative term
    ey=y-k2.*d2+k4.*d4; % ey=sharpened signal
    TwoWarea1b=sum(y(val2ind(x,pos1-FWHM):(val2ind(x,pos1+FWHM)))).*dx;
    TwoWarea2b=sum(y(val2ind(x,pos2-FWHM):(val2ind(x,pos2+FWHM)))).*dx;
    TwoWarea1=sum(ey(val2ind(x,pos1-FWHM):(val2ind(x,pos1+FWHM)))).*dx;
    TwoWarea2=sum(ey(val2ind(x,pos2-FWHM):(val2ind(x,pos2+FWHM)))).*dx;
    
% Remove commenting below to plot animation of signals with different sharpening values
%     figure(1)
%     clf
%     plot(x,y,x,ey)
%     axis([min(x) max(x) 0 2.5])
%     xlabel(['Blue=Original Gaussians, resolution, Rs = ' num2str(Resolution) '    Red = Sharpened    k2 = ' num2str(k2)  ])
%     ylabel('signal amplitude')
%     title('Effect of height ratio on normal and sharpened area measurements of two overlapping Gaussians')
%     
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
        disp(['Regular perpendicular drop failed at amp2 = ' num2str(amp2)])
        PDError1b(trial)=NaN;
        PDError2b(trial)=NaN;
    end
    try
        APresults=autopeaks(x,ey,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,3);PerpDrop=APresults(:,5);
        PDError1(trial)=100*(PerpDrop(1)-TrueArea1)/TrueArea1;
        PDError2(trial)=100*(PerpDrop(2)-TrueArea2)/TrueArea2;
    catch
        disp(['Sharpened perpendicular drop failed at amp2 = ' num2str(amp2)])
        PDError1(trial)=NaN;
        PDError2(trial)=NaN;
    end
    valley(trial)=min(ey(val2ind(x,pos1):val2ind(x,pos2)));
end

% Plot signals with optimum sharpening for minimum errors
figure(1)
clf
plot(x,y,x,ey)
axis([11 25 -0.1 2]);
xlabel(['Blue=Original Gaussians, Resolution, Rs = ' num2str(Rs) '     factor = ' num2str(factor) '    Red: Sharpened peaks' ])
ylabel('signal amplitude')
title('Effect or peak height ratio on errors of perpendicular drop area measurements')
% trial=val2ind(PDError2,min(abs(PDError2)));
text(pos1,amp1,['    Normal error = ' num2str(PDError1b(trial))  '%'] )
text(pos1,max(ey),['    Sharpened error = ' num2str(PDError1(trial))  '%'] )
text(pos2,amp2,['    Normal error = ' num2str(PDError2b(trial)) '%'] )
text(pos2,ey(val2ind(x,pos2)),['    Sharpened error = ' num2str(PDError2(trial)) '%'] )

figure(2)
% AverageError=(abs(PDError1)+abs(PDError2))/2;
% AverageErrorb=(abs(PDError1b)+abs(PDError2b))/2; % average absolute error
plot(amp2vector,abs(PDError2b),'b',amp2vector,abs(PDError2),'r')
axis([0 2 0 5])
title(['Blue: Original peak.   Red: Sharpened peak.   Rs = ' num2str(Rs) '     Factor = ' num2str(factor) ])
ylabel('Absolute percent peak 2 area error')
xlabel('Relative height of peak 2')
disp(['Sharpening factor = ' num2str(factor) ] )
disp(['Calculated K2 = ' num2str((1.7*(FWHM).^2)/factor) ] )
disp(['Calculated K4 = ' num2str(1.7*((FWHM).^4)/(28*factor)) ] )