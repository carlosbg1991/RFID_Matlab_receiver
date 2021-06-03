format compact
format short g
x=14:.01:22;
y=exp(-(x-18).^2);
a=-1300; % Adjust this constant for best results
b=1300000; % Adjust this constant for best results
c=-700000000; % Adjust this constant for best results
d2n=a.*deriv2(y); % Second derivative term
d4=b.*deriv4(y); % Fourth derivative term
d6=c.*deriv4(deriv2(y)); % Sixth derivative term
ey=y+d2n+d4+d6;
figure(1)
plot(x,y,x,d2n,x,d4,x,d6,x,ey)
xlabel(['Derivative weighting factors:    2nd =  ' num2str(a) '    4th = ' num2str(b)  '    6th = ' num2str(c)] )
ylabel('signal or derivative amplitude')
title('Blue=Original Gaussian. Green=Negative 2nd derivative.  Red=4th derivative.  Magenta=6th derivative.  Cyan=Sum of all')
FWHMy=halfwidth(x,y);
FWHMey=halfwidth(x,ey);
yArea=trapz(x,y);
eyArea=trapz(x,ey);
text(18,1,['   FWHM =  ' num2str(FWHMy) '    Area = ' num2str(yArea)] )
text(18,max(ey),['   FWHM =  ' num2str(FWHMey) '    Area = ' num2str(eyArea)] )
PercentChangeInPeakWidth=100.*(FWHMy-FWHMey)/FWHMy;
PercentChangeInPeakArea=100.*(yArea-eyArea)/eyArea;
bl=ey(1:250); % measure baseline flatness
PercentBaselineMax=100.*(max(bl)./max(ey));
disp(' ')
disp('           a        b             c       %ChangeWidth   %ChangeArea  % Baseline')
disp([a b c PercentChangeInPeakWidth PercentChangeInPeakArea PercentBaselineMax])
