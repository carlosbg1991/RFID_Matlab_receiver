% Script that shows how to optimize first derivative symmetrization on two
% overlapping exponentially broadened Gaussians. Plots and compares
% original and sharpened peaks (Figure 1). Tries first-derivative weighting
% factors from +20% to -20% of the correct tau value in line 15 and finds
% the factor that gives the lowest error for peak 2 (Figure 2). Change the
% resolution by changing either the peak positions in lines 20 and 21 or
% the peak width in line 15. Change peak 2 height in line 19. You must have
% derivxy.m, autopeaks.m, val2ind.m and halfwidth.m in the path.
format compact
format short g
clf
scale=1; % x-axis scaling factor
dx=0.01*scale;  % Sampling interval
x=11*scale:dx:25*scale;  % time axis, CHANGE as desired.
w=1*scale; % FWHM (full width at half maximum) of isolated single peak
sigma=w./2.355; % Sigma (standard deviation) of underlying Gaussian
tau=.5; % Correct tau value, <<< CHANGE as desired.
amp1=.99; % Amplitude of peak 1.
amp2=.3;% Amplitude of peak 2, <<< CHANGE as desired.
pos1=13; % Position (retention time) of peak 1.
pos2=15; % Position (retention time) of peak 2. CHANGE as desired.
y1=amp1.*expgaussian(x,pos1*scale,w,-tau./dx);
FWHM=halfwidth(x,y1,pos1);
y2=amp2.*expgaussian(x,pos2*scale,w,-tau./dx);
y=y1+y2;
TrueArea1=trapz(x,y1);
TrueArea2=trapz(x,y2);
Resolution=(pos2-pos1)/(4*FWHM/2.355); % Estimate of 4*sigma resolution for broadened peak
% disp('        factor     FWHMbefore   FWHMafter    %PDError1   %PDError2 ')

% Calculations
for trial=1:80
    factor(trial)=.5*tau+trial*.01.*tau;
    d=derivxy(x,y); % First derivative term
    sy=y+factor(trial).*d; % sy=trial sharpened signal for this factor
    TwoWarea1b=sum(y(val2ind(x,pos1-w):(val2ind(x,pos1+w)))).*dx;
    TwoWarea2b=sum(y(val2ind(x,pos2-w):(val2ind(x,pos2+w)))).*dx;
    TwoWarea1=sum(sy(val2ind(x,pos1-w):(val2ind(x,pos1+w)))).*dx;
    TwoWarea2=sum(sy(val2ind(x,pos2-w):(val2ind(x,pos2+w)))).*dx;
    
    xlabel(['Blue=Original signal, Resolution, Rs = ' num2str(Resolution) '     Red=Sharpened signal, factor '  num2str(factor(trial))  ])
    ylabel('signal amplitude')
    title('Effect of sharpening on area measurements of overlapping exp. broadened Gaussians')
    FWHMy1=halfwidth(x,y,pos1);
    FWHMey1=halfwidth(x,sy,pos1);
    FWHMy2=halfwidth(x,y,pos2);
    FWHMey2=halfwidth(x,sy,pos2);
    yArea=trapz(x,y);
    eyArea=trapz(x,sy);
    PercentChangeInPeak1Width=100.*(FWHMy1-FWHMey1)/FWHMy1;
    PercentAreaError1b=100*(TrueArea1-TwoWarea1b)/TrueArea1;
    PercentAreaError2b=100*(TrueArea2-TwoWarea2b)/TrueArea2;
    PercentAreaError1=100*(TrueArea1-TwoWarea1)/TrueArea1;
    PercentAreaError2=100*(TrueArea2-TwoWarea2)/TrueArea2;
    
    SlopeThreshold=.0000000001;
    AmpThreshold=.01;
    smoothwidth=3;
    peakgroup=5;
    
    try
        APresultsb=autopeaks(x,y',SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,3);PerpDropb=APresultsb(:,5);
        PDError1b(trial)=100*(PerpDropb(1)-TrueArea1)/TrueArea1;
        PDError2b(trial)=100*(PerpDropb(2)-TrueArea2)/TrueArea2;
    catch
        PDError1b(trial)=NaN;
        PDError2b(trial)=NaN;
    end
    
    try
        APresults=autopeaks(x,sy',SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,3);PerpDrop=APresults(:,5);
        PDError1(trial)=100*(PerpDrop(1)-TrueArea1)/TrueArea1;
        PDError2(trial)=100*(PerpDrop(2)-TrueArea2)/TrueArea2;
    catch
        PDError1(trial)=NaN;
        PDError2(trial)=NaN;
    end
    
   %  disp([factor(trial) FWHMy1 FWHMey1 PDError1(trial) PDError2(trial)])
    figure(1)
    plot(x,y,x,sy)
    text(pos1,max(y1),['    Area error = ' num2str(PDError1b(trial))  '%'] )
    text(pos1,max(sy),['    Area error = ' num2str(PDError1(trial))  '%'] )
    text(pos2,max(y2),['    Area error = ' num2str(PDError2b(trial)) '%'] )
    text(pos2,sy(val2ind(x,pos2)),['    Area error = ' num2str(PDError2(trial)) '%'] )
    xlabel(['Blue=Original signal, Resolution, Rs = ' num2str(Resolution) '     Red=Sharpened signal, factor '  num2str(factor(trial))  ])
    ylabel('signal amplitude')
    title('Effect of symmetrization on area measurements of overlapping exp. broadened Gaussians')
    drawnow
end
OptIndex1=val2ind(PDError1,min(abs(PDError1)));
OptIndex2=val2ind(PDError2,min(abs(PDError2)));
OptFactor1= factor(OptIndex1);
OptFactor2= factor(OptIndex2);

figure(2)
plot(factor,abs(PDError1),factor,abs(PDError2),factor,abs(PDError1b),factor,abs(PDError2b))
title(['Blue = Symm. peak 1     Red = Symm. peak 2     Resolution, Rs = ' num2str(Resolution) '     Tau = '  num2str(tau)  ])
ylabel('Absolute percent area error')
xlabel('First derivative weighting factor')

figure(1)
plot(x,y,x,sy)
text(pos1,max(y1),['    Area error = ' num2str(PDError1b((OptIndex2)))  '%'] )
text(pos1,max(sy),['    Sharpened Area error = ' num2str(PDError1((OptIndex2)))  '%'] )
text(pos2,max(y2),['    Area error = ' num2str(PDError2b((OptIndex2))) '%'] )
sy=y+OptFactor2.*d; % sy=sharpened signal
text(pos2,sy(val2ind(x,pos2)),['     Sharpened Area error = ' num2str(PDError2((OptIndex2))) '%'] )
xlabel(['Blue: Original signal, Resolution, Rs = ' num2str(Resolution) ', Tau = ' num2str(tau)  '.    Red: Sharpened signal, factor '  num2str(OptFactor2) ])
ylabel('signal amplitude')
title('Optimum symmetrization for area measurements of overlapping exp. broadened Gaussians')

disp(' ')