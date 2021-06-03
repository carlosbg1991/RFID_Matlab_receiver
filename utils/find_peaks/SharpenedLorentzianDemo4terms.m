format compact
format short g
x=10:.01:26;
y=lorentzian(x,18,2);
d2n=-9500.*deriv2(y); % Adjust this constant for best results
d4=1300000.*deriv4(y); % Adjust this constant for best results
d6=-1000000000.*deriv4(deriv2(y)); % Adjust this constant for best results
ey=y+d2n+d4+d6;
figure(1)
clf
plot(x,y,x,d2n,x,d4,x,d6,x,ey);
axis([10.2 25.8 -.5 4])
xlabel('x')
ylabel('signal or derivative amplitude')
title('Blue=Original Gaussian.  Green=Negative 2nd derivative.  Red=4th derivative.  Cyan=6th derivative  Magenta=Sum of all')
FWHMy=halfwidth(x,y);
FWHMey=halfwidth(x,ey);
yArea=trapz(x,y);
eyArea=trapz(x,ey);
text(18,1,['   FWHM =  ' num2str(FWHMy) '    Area = ' num2str(yArea)] )
text(18,max(ey),['   FWHM =  ' num2str(FWHMey) '    Area = ' num2str(eyArea)] )
PercentChangeInPeakWidth=100.*(FWHMy-FWHMey)/FWHMy
PercentChangeInPeakArea=100.*(yArea-eyArea)/eyArea


