% Demonstrates 2nd plus 4th derivative sharpening on Lorentzian. Plots and
% compares original peak, 2nd and 4th derivatives, and sharpened peak.You
% can adjust "k2" and "k4" in lines 10 and 11. 
% Must have lorentzian.m, derivxy.m and halfwidth.m in the path.
format compact
format short g
scale=1; % x-axis scaling factor
dx=0.1*scale;
x=9*scale:dx:27*scale;
w=3*scale;
y=lorentzian(x,18*scale,w);
% Estimated values of weighting factors k2 and k4 based on peak width
k2=((w).^2)/4; % <<<< You can adjust this 2nd derivative constant for best results 
k4=((w).^4)/600; % <<<< You can adjust this 4th derivative constant for best results 
d=derivxy(x,y); % First derivative term
d2=derivxy(x,d);% Second derivative term
d3=derivxy(x,d2);% Third derivative term
d4=derivxy(x,d3); % Fourth derivative term
ey=y-k2.*d2+k4.*d4; % ey=sharpened signal
figure(1)
clf
plot(x,y,x,ey,x,-k2.*d2,x,k4.*d4)
grid
% axis([10.2 25.8 -.5 2])
xlabel(['Derivative weighting factors:    2nd =  ' num2str(k2) '    4th = ' num2str(k4)  ] )
ylabel('signal or weighted derivative amplitudes')
title('Blue=Original lorentzian.   Yellow=Negative 2nd derivative.   Purple=4th derivative.   Red=Sum of all')
FWHMy=halfwidth(x,y);
FWHMey=halfwidth(x,ey);
yArea=trapz(x,y);
eyArea=trapz(x,ey);
text(18,1,['   FWHM =  ' num2str(FWHMy) '    Area = ' num2str(yArea)] )
text(18,max(ey),['   FWHM =  ' num2str(FWHMey) '    Area = ' num2str(eyArea)] )
PercentChangeInPeakWidth=100.*(FWHMy-FWHMey)/FWHMy;
PercentChangeInPeakArea=100.*(yArea-eyArea)/eyArea;
disp(' ')
%disp('          k2        k4        %ChangeWidth  %ChangeArea  ')
% disp([k2 k4 PercentChangeInPeakWidth PercentChangeInPeakArea ])

OriginalTotalArea=yArea;
SharpenedTotalArea=eyArea;
figure(2)
[FitResults,GOF]=peakfit([x;ey],0,0,1,1);
disp('        k2           k4          FWHM      Gaussian error')
disp([k2 k4 FWHMey GOF(1)])