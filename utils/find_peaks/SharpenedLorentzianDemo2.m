% Demonstrates 2nd plus 4th derivative sharpening followed by power method.
% Plots and compares original peak, even-derivative sharpened peak, and the
% effect of additional power method. You can adjust "k2" and "k4" in lines
% 18 and 19. Shows that the fraction of total area in the central maximum
% is about 99% Must have lorentzian.m, derivxy.m, peakfit.m, and halfwidth.m
% in the path.
format compact
format short g
scale=1; % x-axis scaling factor
dx=0.01*scale;
x=12*scale:dx:24*scale;
lx=length(x);
w=1.2*scale;
Noise=0.0;
y=lorentzian(x,18*scale,w)+Noise.*randn(size(x));
SmoothWidth=51;
y=fastsmooth(y,SmoothWidth,3,1);
% Estimated values of weighting factors k2 and k4 based on peak width
k2=((w).^2)/5; % <<<< You can adjust this 2nd derivative constant for best results 
k4=((w).^4)/100; % <<<< You can adjust this 4th derivative constant for best results 
% k4=0; % Second derivative sharpening only in this example
d=derivxy(x,y); % First derivative term
d2=derivxy(x,d);% Second derivative term
d3=derivxy(x,d2);% Third derivative term
d4=derivxy(x,d3); % Fourth derivative term
ey=(y-k2.*d2+k4.*d4); % ey=sharpened signal
eyy=ey.*(y./max(y)); % multiplied by y/max(y)
% Clean up edges
x=x(SmoothWidth:lx-SmoothWidth);
y=y(SmoothWidth:lx-SmoothWidth);
ey=ey(SmoothWidth:lx-SmoothWidth);
eyy=eyy(SmoothWidth:lx-SmoothWidth);
ey2=ey.^2;
figure(1)
clf
plot(x,y,x,ey,x,eyy,'g',x,ey2,'c')
grid
% axis([10.2 25.8 -.5 2])
xlabel(['Derivative weighting factors:    2nd =  ' num2str(k2) '    4th = ' num2str(k4)  ] )
ylabel('signal or weighted derivative amplitudes')
title('Blue=Original, y.  Red=Even deriv sharp, ey.  Green=ey x y/max(y), eyy.  Cyan=ey^2.')
FWHMy=halfwidth(x,y);
FWHMey=halfwidth(x,ey);
FWHMeyy=halfwidth(x,eyy);
FWHMey2=halfwidth(x,ey2);

yArea=trapz(x,y);
eyArea=trapz(x,ey);
eyyArea=trapz(x,eyy);
ey2Area=trapz(x,ey2);

text(18,max(y),['   FWHM =  ' num2str(FWHMy) '    Area = ' num2str(yArea)] )
text(18,max(ey),['   FWHM =  ' num2str(FWHMey) '    Area = ' num2str(eyArea)] )
text(18,max(eyy),['   FWHM =  ' num2str(FWHMeyy) '    Area = ' num2str(eyyArea)] )

PCIPeakWidth=100.*(FWHMy-FWHMey)/FWHMy;
PCIPeakArea=100.*(yArea-eyArea)/eyArea;
PCIPeakWidth2=100.*(FWHMy-FWHMeyy)/FWHMy;
PCIPeakArea2=100.*(yArea-eyyArea)/eyyArea;
PCIPeakWidth3=100.*(FWHMy-FWHMey2)/FWHMy;
PCIPeakArea3=100.*(yArea-ey2Area)/eyyArea;
disp(' ')
%disp('          k2          k4      %ChangeWidth  %ChangeArea  ')
% disp([k2 k4 PCIPeakWidth PCIPeakArea ])

OriginalTotalArea=yArea;
SharpenedTotalArea=eyArea;

figure(2);[FitResults,GOF]=peakfit([x;ey],0,0,1,1);
subplot(2,1,1)
title('Gaussian fit to even-derivative sharpened, ey')

figure(3);[FitResults2,GOF2]=peakfit([x;eyy],0,0,1,1);
subplot(2,1,1)
title('Gaussian fit to even-derivative sharpened * y, eyy')

figure(4);[FitResults3,GOF3]=peakfit([x;ey.^2],0,0,1,1);
subplot(2,1,1)
title('Gaussian fit to power-raised even-derivative sharpened, ey^2')

disp('Sharpened:   % Width Change  % Area Change   Gaussian R2')
disp(['                 ' num2str([PCIPeakWidth PCIPeakArea GOF(2)] )] )
disp('Sharpened*y: % Width Change  % Area Change   Gaussian R2')
disp(['                 ' num2str([PCIPeakWidth2 PCIPeakArea2 GOF2(2)]) ])
disp('Sharpened^2: % Width Change  % Area Change   Gaussian R2')
disp(['                 ' num2str([PCIPeakWidth3 PCIPeakArea3 GOF3(2)]) ])