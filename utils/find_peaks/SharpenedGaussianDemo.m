% Demonstrates 2nd plus 4th derivative sharpening on Gaussian.
% Plots and compares original peak, 2nd and 4th derivatives, and
% sharpened peak.You can adjust "k2" and "k4" in lines 10 and 11.
% Shows that the fraction of total area in the central maximum is about 99%
% Must have gaussian.m, derivxy.m, peakfit.m, and halfwidth.m in the path.
format compact
format short g
scale=1; % x-axis scaling factor
dx=0.01*scale;
x=12*scale:dx:24*scale;
w=2*scale;
y=gaussian(x,18*scale,w);
% Estimated values of weighting factors k2 and k4 based on peak width
k2=((w).^2)/32; % <<<< You can adjust this 2nd derivative constant for best results 
k4=((w).^4)/900; % <<<< You can adjust this 4th derivative constant for best results 
k4=0;
d=derivxy(x,y); % First derivative term
d2=derivxy(x,d);% Second derivative term
d3=derivxy(x,d2);% Third derivative term
d4=derivxy(x,d3); % Fourth derivative term
ey=(y-k2.*d2+k4.*d4).*(y./max(y)).*sqrt(2); % ey=sharpened signal
figure(1)
clf
plot(x,y,x,ey,x,-k2.*d2,x,k4.*d4)
% axis([10.2 25.8 -.5 2])
xlabel(['Derivative weighting factors:    2nd =  ' num2str(k2) '    4th = ' num2str(k4)  ] )
ylabel('signal or weighted derivative amplitudes')
title('Blue=Original Gaussian.   Yellow=Negative 2nd derivative.   Purple=4th derivative.   Red=Sum of all')
FWHMy=halfwidth(x,y);
FWHMey=halfwidth(x,ey);
yArea=trapz(x,y);
eyArea=trapz(x,ey);
text(18,1,['   FWHM =  ' num2str(FWHMy) '    Area = ' num2str(yArea)] )
text(18,max(ey),['   FWHM =  ' num2str(FWHMey) '    Area = ' num2str(eyArea)] )
PercentChangeInPeakWidth=100.*(FWHMy-FWHMey)/FWHMy;
PercentChangeInPeakArea=100.*(yArea-eyArea)/eyArea;
disp(' ')
%disp('          k2          k4      %ChangeWidth  %ChangeArea  ')
% disp([k2 k4 PercentChangeInPeakWidth PercentChangeInPeakArea ])

OriginalTotalArea=yArea;
SharpenedTotalArea=eyArea;
figure(2);[FitResults,GOF]=peakfit([x;ey],0,0,1,1);

[FitResults,GOF]=peakfit([x;ey],0,0,1,1,0,0,0,0,0,0);
disp('        k2           k4       Width Change  AreaChange   Gaussian %error')
disp([k2 k4 PercentChangeInPeakWidth PercentChangeInPeakArea GOF(1)])