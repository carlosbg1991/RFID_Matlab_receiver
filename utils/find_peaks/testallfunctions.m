% Testallfunctions test most functions and scripts listed on  
% http://terpconnect.umd.edu/~toh/spectrum/functions.html
% Approimately 200 scripts and functions tested as of February 2018.
clf
disp('Tests most functions and scripts. Takes about 4 to 40 minutes to')
disp('complete in Matlab, depending on the speed of your computer.')
disp('Leaves a text file called diary with list of all items tested.')
diary off

format short g
format compact

clear;diary on;TestStartTime=cputime
save TestStartTime
diary off
fn=0;
save fn
clear;load fn;load fn;fn=fn+1;save fn;save fn;diary on;disp([num2str(fn) ': alphafunction(x,pos,spoint) ']);diary off
x=0:.1:25;y=alphafunction(x,2,8);plot(x,y)

clear;load fn;load fn;fn=fn+1;save fn;save fn;diary on;disp([num2str(fn) ':bimodal([1:1000],.1,-.2,.2) ']);diary off
hist(bimodal(1:1000,.1,-.2,.2),30)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':g = BiGaussian(x,pos,wid,m)  ']);diary off
x=0:.1:25;y=BiGaussian(x,10,8,50);plot(x,y)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':g = BiLorentzian(x,pos,wid,m) ']);diary off
x=0:.1:25;y=BiLorentzian(x,10,8,50);plot(x,y)

clear;load fn;fn=fn+1;save fn;diary on;disp([ num2str(fn) ': bmd=bimodal(x,std,a,b) ']);diary off
hist(bimodal(1:1000,.1,-.2,.2),30);

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':BlackbodyDataFit ']);diary off
BlackbodyDataFit

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':bluenoise(NumPoints)  ']);diary off
plot(bluenoise(1000))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':[Height,Position,Width,BootResults]=bootgaussfit(x,y,plots)  ']);diary off
x=50:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
[Height,Position,Width,BootResults]=bootgaussfit(x,y,1);

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':BootstrapIterativeFit(100,100,100,20,10,100)  ']);diary off
clf;BootstrapIterativeFit(100,100,100,20,10,100)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':BootstrapIterativeFit2(100,100,100,50,200,100,100,10,100);  ']);diary off
clf;BootstrapIterativeFit2(100,100,100,50,200,100,100,10,100);

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':BWF (Breit-Wigner-Fano)   ']);diary off
clf;x=0:.1:25;y=BWF(x,10,8,5);plot(x,y)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':CaseStudyC  ']);diary off
CaseStudyC
pause(.1)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':CentralLimitDemo  ']);diary off
CentralLimitDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':clippedgaussian  ']);diary off
x=0:.1:25;y=clippedgaussian(x,10,8,5);plot(x,y);diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':clippedlorentzian  ']);diary off
x=0:.1:25;y=clippedlorentzian(x,10,8,.8);plot(x,y);diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':cls ']);diary off
clsdemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':CLSvsINLS  ']);diary off
CLSvsINLS

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':sy=condense(y,n)  ']);diary off
condense([1 2 3 4 5 6 7 8 9 10 11 12],3)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':sm=condensem(M,n)  ']);diary off
condensem([1 2 3 4 ; 4 5 6 7],2)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':CurvefitNoiseColorTest  ']);diary off
CurvefitNoiseColorTest

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DeconvDemo  ']);diary off
DeconvDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DeconvDemo2  ']);diary off
DeconvDemo2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':deconvolutionexample  ']);diary off
deconvolutionexample

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Demo3peaks  ']);diary off
Demo3peaks

clf
clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoFindPeak ']);diary off
DemoFindPeak;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoFindPeaksb ']);diary off
DemoFindPeaksb;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoFindPeaksb3 ']);diary off
DemoFindPeaksb3;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoFindPeakSNR ']);diary off
DemoFindPeakSNR;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoFindsquare ']);diary off
if IsOctave,pkg load signal,end
DemoFindsquare;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Demofitblackbody ']);diary off
Demofitblackbody;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Demofitgauss (Single Gaussian) ']);diary off
Demofitgauss;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Demofitgauss2 (two overlapping Gaussians) ']);diary off
Demofitgauss2;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Demofitgaussb ']);diary off
Demofitgaussb;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Demofitlorentzianb ']);diary off
Demofitlorentzianb;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Demofitmultiple ']);diary off
Demofitmultiple;

clear
clf
clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':demoipeak ']);diary off
demoipeak;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Demoipf ']);diary off
if IsOctave,
else
  Demoipf
end

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':demoisignal ']);diary off
demoisignal;


clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':deriv ']);diary off
a=[1 1 1 2 3 5 6 8 9 10 10 10 10 10];
deriv(a)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':deriv1 ']);diary off
a=[1 1 1 2 3 5 6 8 9 10 10 10 10 10];
deriv1(a)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':deriv2 ']);diary off
a=[1 1 1 2 3 5 6 8 9 10 10 10 10 10];deriv2(a)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':deriv3 ']);diary off
a=[1 1 1 2 3 5 6 8 9 10 10 10 10 10];deriv3(a)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':deriv4 ']);diary off
a=[1 1 1 2 3 5 6 8 9 10 10 10 10 10];deriv4(a)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DerivativeDemo ']);diary off
DerivativeDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DerivativeShapeDemo ']);diary off
DerivativeShapeDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':derivdemo1 ']);diary off
derivdemo1

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':derivdemo2  ']);diary off
derivdemo2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': derivdemo3 ']);diary off
 derivdemo3

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':derivdemo4  ']);diary off
derivdemo4

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':derivxy ']);diary off
x = [1 2 4 7 10 14];
y = 2*x;
derivxy(x,y)

clf
clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':downsigmoid  ']);diary off
plot(downsigmoid(1:100,50,50))

clf
clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':enhance(signal,factor1,factor2,SmoothWidth) ']);diary off
   x=0:.01:18;
   y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2);
   Enhancedsignal=enhance(y,1000,1000000,3);
   plot(x,y,x,Enhancedsignal,'r ');diary off
   xlabel('Blue: original signal.    Red: enhanced signal. ');diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ExpBroaden  ']);diary off
 ExpBroaden([0 0 0 1 1 2 4 6 4 3 2 1 1 0 0 0 ]',-2)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ExpBroadenDemo  ']);diary off
ExpBroadenDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':expgaussian  ']);diary off
x=0:.01:18;
plot(x,gaussian(x,6,2),x,expgaussian(x,6,2,-30),'r ');diary off
xlabel('Blue: Gaussian.    Red: exponentially broadened  Gaussian. ');diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':explorentzian  ']);diary off
x=0:.01:18;
plot(x,lorentzian(x,6,2),x,explorentzian(x,6,2,-30),'r ');diary off
xlabel('Blue: Lorentzian.    Red: exponentially broadened Lorentzian. ');diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':exppulse  ']);diary off
x=1:100;
plot(x,exppulse(x,10,15))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':fastsmooth  ']);diary off
x=1:100;
y=randn(size(x)); 
plot(x,y,x,fastsmooth(y,5,3,1),'r ');diary off
xlabel('Blue: white noise.    Red: smoothed white noise. ');diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksG ']);diary off
findpeaksG(0:.01:2,humps(0:.01:2),0,-1,5,5)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksb  ']);diary off
x=1:.2:100;
y=modelpeaks(x,4,2,[10 2 2 2],[-10 20 50 80],[50 3 3 3])+.03*randn(size(x));
disp('          Peak      Position      Height      Width        Area         % error      R2 ');
NoBackgroundSubtraction=findpeaksplot(x,y,.00005,.5,15,15,3)
QuadraticBackgroundSubtraction=findpeaksb(x,y,.00005,.5,15,15,3,150,2,0,2,2)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksb3  ']);diary off
x=1:.2:100;
y=modelpeaks(x,4,2,[1 2 3 4],[45 50 55 60],[3 3 3 3])+.05*randn(size(x));
findpeaksresults=findpeaksplot(x,y,.00005,.5,9,12,3);
findpeaksb3results=findpeaksb3(x,y,.00005,.5,9,12,3,2,0,10,0,0);
disp('          Peak      Position      Height      Width        Area        % error ');
findpeaksError=mean(findpeaksresults(:,3)-[1 2 3 4] );
findpeaksb3Error=mean(findpeaksb3results(:,3)-[1 2 3 4] );

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksE ']);diary off
x=[0:.01:50];
y=cos(x);
P=findpeaksE(x,y,0,-1,5,5)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksfit  ']);diary off
x=1:.2:100;
y=lorentzian(x,10+30.*rand(),10)+lorentzian(x,50,10)+lorentzian(x,60+30.*rand(),10)+.05.*randn(size(x));
[findpeaksr,peakfitr]=findpeaksfit(x,y,8e-005,0.5,17,31,3,2,0,1,0,0,1);
FindPeaksError=mean(findpeaksr(:,3)-[1 1 1] );diary off
FindPeaksFitError=mean(peakfitr(:,3)-[1 1 1] );diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksGSS  ']);diary off
  x=1:.2:100;
  y=gaussian(x,20,10)+gaussian(x,50,10)+gaussian(x,80,10)+.1.*randn(size(x));
  disp('           Peak      Position     Height      Width        Area       Start        End ');
  findpeaksGSS(x,y,0.0004,0.3,17,21,3)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksL  ']);diary off
    x=1:.2:100;
    y=lorentzian(x,20,10)+lorentzian(x,50,10)+lorentzian(x,80,10)+.01.*randn(size(x));
    findpeaksL(x,y,0.0004,0.3,17,17,3)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksLSS  ']);diary off
x=1:.2:100;
y=lorentzian(x,20,5)+lorentzian(x,50,5)+lorentzian(x,80,5)+.1.*randn(size(x));
findpeaksLSS(x,y,0.0004,0.3,17,21,3)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksnr  ']);diary off
findpeaksnr(0:.01:2,humps(0:.01:2),0,-1,5,5)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksplot  ']);diary off
findpeaksplot(0:.01:2,humps(0:.01:2),0,-1,5,5)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksplotL  ']);diary off
findpeaksplotL(0:.01:2,humps(0:.01:2),0,-1,5,5)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksT  ']);diary off
x=[-10:.2:10];
y=exp(-(x+5*randn()).^2)+.01.*randn(size(x));
findpeaksT(x,y,0.003,0.5,7,9,3)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findpeaksTplot  ']);diary off
x=-10:.2:10;
y=exp(-(x+5*randn()).^2)+.01.*randn(size(x));
findpeaksTplot(x,y,0.003,0.5,7,9,3)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findsquarepulse  ']);diary off
t=1:100;
y=sign(sin(t));plot(t,y);
S=findsquarepulse(t,y,0)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':findvalleys  ']);diary off
x=[0:.01:50];y=cos(x);P=findvalleys(x,y,0,-1,5,5)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':fitblackbody  ']);diary off
clf;BlackbodyDataFit

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':fitgauss2  ']);diary off
Demofitgauss

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':fitgauss2b  ']);diary off
Demofitgaussb

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':fitlorentzianb  ']);diary off
Demofitlorentzianb

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':fitM  ']);diary off
options = optimset('TolX',0.000001);
absorbance=fminsearch(@(lambda)(fitM(lambda,[0.56529 0.38696 0.56529 0.73496]',[0.2 1 0.2 0.058824]',[1 0.5 0.0625 0.5]',.01)),1)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':fitmultiple  ']);diary off
Demofitmultiple

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':FouFilter  ']);diary off
clf;subplot(211);y=randn(size(1:1000));plot(y);
title('White noise ');
subplot(212);plot(FouFilter(y,1,19,2,2,0))
title('narrow band of frequencies ');


clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':gaussfit  ']);diary off
[Height, Position, Width]=gaussfit([1 2 3],[1 2 1]) 

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':GL (Gaussian/Lorentzian blend)  ']);diary off
clf;x=1:100;plot(GL(x,50,10,50))
title('50% Gaussian, 50% Lorentzian ');diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':idpeaks  ']);diary off
warning offfindpeaksb3
load spectrum  % Cu atomic emission spectrum
load DataTableCu  % contains Positions,Names for Cu atomic lines
IdentifiedPeaks=idpeaks(Cu,0.01,.001,5,5,.01,Positions,Names);
for number=2:10;
   disp(['wavelength=' num2str([cell2mat(IdentifiedPeaks(number,1))]) '      element=' num2str([cell2mat(IdentifiedPeaks(number,2))]) ] )
end

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ifilter  ']);diary off
x=0:1000; 
y=sin(x/10)+sin(x/5)+sin(x);
ifilter(x,y);

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ipeak  ']);diary off
x=[0:.1:100];
y=(sin(x)).^2+.1*randn(size(x));
ipeak(x,y);

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ipeakdemo  ']);diary off
ipeakdemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ipeakdemo1  ']);diary off
ipeakdemo1

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ipeakdemo2  ']);diary off
ipeakdemo2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ipeakdemo3  ']);diary off
ipeakdemo3

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ipeakdemo4  ']);diary off
ipeakdemo4

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ipeakdemo5  ']);diary off
ipeakdemo5

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ipeakdemo6  ']);diary off
ipeakdemo6

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':isignaldemo  ']);diary off
isignaldemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':IsOctave  ']);diary off
IsOctave

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':LeastSquaresMatlab ']);diary off
x=1:10;y=[1 3 4 2 5 6 4 9 8 9];
clf
LeastSquaresMatlab

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':LinearFiMC  ']);diary off
LinearFiMC

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':lognormal  ']);diary off
x=0:.1:10;plot(x,lognormal(x,5,20))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':lorentzfit  ']);diary off
x=50:150;y=100.*lorentzian(x,100,100)+10.*randn(size(x));
[Height,Position,Width]=lorentzfit(x,y)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':lorentzian ']);diary off
x=0:.1:10;plot(x,lorentzian(x,5,2))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':medianfilter  ']);diary off
x=-5:.001:5;y=exp(-(x).^2);
for n=1:1000:10000,y(n)=randn(1)+y(n);end
subplot(1,2,1);plot(x,y);title('Before ');diary off
subplot(1,2,2);plot(x,medianfilter(y,1));title('After ');diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':modelpeaks ']);diary off
clf
x=1:1000;
model=modelpeaks(x,3,1,[1 2 3],[200 500 700],[50 50 50],0,0);
plot(model)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':modelpeaks2 ']);diary off
clf
t=-10:.1:50;
plot(t,modelpeaks2(t,[1 2],[2 1],[10 30],[5 5],[0 0]))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':MonteCarloDemo  ']);diary off
MonteCarloDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ngaussian  ']);diary off
x=0:.1:10;plot(x,ngaussian(x,5,2,2))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':NoiseColorTest  ']);diary off
clf;NoiseColorTest

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':noisetest  ']);diary off
clf;noisetest

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':NumPeaksTest  ']);diary off
clf;NumPeaksTest

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':peakfit  ']);diary off
x=0:.1:10;
y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
peakfit([x' y'],0,0,2)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PeakFitDemo11Gauss  ']);diary off
PeakFitDemo11Gauss

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PeakFitDemo11Lor  ']);diary off
PeakFitDemo11Lor

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':peakfitdemob  ']);diary off
peakfitdemob

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PeakfitTimeTest  ']);diary off
PeakfitTimeTest

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PeakfitTimeTest2  ']);diary off
PeakfitTimeTest2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PeakfitTimeTest3  ']);diary off
PeakfitTimeTest3

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':peakfitVSfindpeaks  ']);diary off
peakfitVSfindpeaks

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': peakfunction ']);diary off
x=1:100;plot(x,peakfunction(4,x,50,30,3))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':peakstats  ']);diary off
clf;
x=[0:.1:100];y=5+5.*cos(x)+randn(size(x));peakstats(x,y,0,-1,11,19,3,1);

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':pearson ']);diary off
clf;
x=1:100;plot(x,pearson(x,50,5,3))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':pinknoise  ']);diary off
clf;
plot(pinknoise(1000));

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':plotfita  ']);diary off
clf;
x=1:100;
y=1+x+5.*randn(size(x));
[coef, RSquared,BootResults]=plotfita(x,y,2);

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PlotFrequencySpectrum  ']);diary off
clf;
x=[0:.01:2*pi]';
y=sin(200*x)+randn(size(x));
subplot(2,1,1);
plot(x,y);
subplot(2,1,2);
PowerSpectrum=PlotFrequencySpectrum(x,y,1,0,1);
title('Power spectrum of very noisy sine wave (SNR=0.7) ');diary off


clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':plotgaussfit  ']);diary off
[Height, Position, Width]=plotgaussfit([0 1 2 3 4],[.1 1 2 1 .1])
title('[Height, Position, Width]=plotgaussfit.m')

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':plotit  ']);diary off
clf;x=[1 2 3 4 5];y=[0 2 3 3 5];plotit(x,y,3,'or','-g ');diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PlotLogPowerSpectrum  ']);diary off
clf;x=[0:.01:2*pi]';PlotLogPowerSpectrum(sin(200*x)+randn(size(x)));

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PlotPowerSpectrum  ']);diary off
clf;x=[0:.01:2*pi]';PlotPowerSpectrum(x,sin(200*x)+randn(size(x)));

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':plotting  ']);diary off
clf;plotting

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':plotting2  ']);diary off
clf;plotting2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': ProcessSignal ']);diary off
clf;
x=[0:.01:10]';
y=10./(1+(x/4).^2)+exp(-(x-5).^2)+.1.*randn(size(x));
subplot(2,1,1);plot(x,y);
subplot(2,1,2);plot(x,ProcessSignal(x,y,2,60,3,0,0,0,0,0,0));
title('Use of smoothed second derivative to reduce influence of background on the weak peak at x=5 ');diary off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':propnoise (proportional noise)  ']);diary off
clf;clear
model=gaussian(1:1000,500,200);
model=model+.1.*propnoise(model);
plot(model)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':QuadFitToGaussian  ']);diary off
QuadFitToGaussian

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':RANDtoRANDN  ']);diary off
RANDtoRANDN

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':RegressionDemo  ']);diary off
RegressionDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':rmnan (Removes all NaN and Inf from vector of numerical data ']);diary off
rmnan([1 2 3 4 NaN 6 7 Inf  9])

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':secderivxy  ']);diary off
x=1:100;
y=gaussian(x,50,20);
plot(x,secderivxy(x,y))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':sft (Slow Fourier Transform)  ']);diary off
x=1:100;plot(sft(sin(x)))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ShapeDemo  ']);diary off
clf
ShapeDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':ShapeTest  ']);diary off
clf;ShapeTest

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':upsigmoid  ']);diary off
x=1:100;plot(upsigmoid(x,1,1))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SmoothExperiment  ']);diary off
SmoothExperiment

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SmoothOptimization  ']);diary off
SmoothOptimization

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SmoothVsFit  ']);diary off
help SmoothVsFit
SmoothVsFit

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SmoothVsFitArea  ']);diary off
help SmoothVsFitArea
SmoothVsFitArea

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SpikeDemo1  ']);diary off
SpikeDemo1

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SpikeDemo2  ']);diary off
SpikeDemo2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':sqrtnoise  ']);diary off
model=gaussian([1:1000],500,200);
model=model+.1.*sqrtnoise(model);
plot(model)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':Testfindpeaksb  ']);diary off
Testfindpeaksb

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':testipeak  ']);diary off
testipeak

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':testisignal ']);diary off
testisignal

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':TestLinearFit  ']);diary off
TestLinearFit

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':TestLinearFitHistogram  ']);diary off
TestLinearFitHistogram

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':testnumpeaks  ']);diary off
x=0:.01:10;
y=exp(-(x-4).^2) + exp(-(x-5.3).^2) + exp(-(x-6.5).^2)+.05.*randn(size(x));
testnumpeaks(x,y,1,0,5,5)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':tfit ']);diary off
clf;tfit(10)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':TFit3 ']);diary off
clf;TFit3([3 .1 5])

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':TFitCalCurve  ']);diary off
TFitCalCurve

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':TFitStats  ']);diary off
TFitStats

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':triangle  ']);diary off
clf
x=1:100;plot(triangle(x,50,30))
title('triangle(x,50,30)')

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':triangulationdemo  ']);diary off
triangulationdemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':triangulationExp  ']);diary off
triangulationExp

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':trydatatrans  ']);diary off
x=1:12;y =exp(x);trydatatrans(x,y,1)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':trypoly  ']);diary off
clf;x=1:12;y =cos(x);bar(trypoly(x,y));

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':trypolyplot  ']);diary off
x = 1:12;y=cos(x);trypolyplot(x,y);

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':tsallis  ']);diary off
clf
help tsallis
hold on;
for q=1.1:.1:2.2;
    x=1:100;plot(tsallis(x,50,10,q));
end;
hold off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':tsmooth  ']);diary off
y=[0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];
plot(tsmooth(y,5))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':upsigmoid ']);diary off
x=1:100;plot(upsigmoid(x,50,30))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':val2ind  ']);diary off
x=[1 2 4 3 5 9 6 4 5 3 1];
val2ind(x,6)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':violetnoise ']);diary off
plot(violetnoise(1:100))

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':voigt  ']);diary off
clf
hold on;
for q=1.1:.1:2.2;
    x=1:100;plot(voigt(x,50,5,q));
end;
hold off

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':VoigtFixedAlpha  ']);diary off
clf;VoigtFixedAlpha

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':VoigtVariableAlpha  ']);diary off
clf;VoigtVariableAlpha

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':wheat  ']);diary off
clf;wheat

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':whitenoise  ']);diary off
clf;plot(whitenoise(1:100))
% 
% 
clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':IQrange ']);diary off
a=randn(1000,1); % Create randome numbers with standard deviation of 1.
std(a) % Display standard deviation
a(500)=100; % Add single outlier
std(a) % Shows that standard deviation is larger with single outlier.
IQrange(a)  % Shows that IQrange is not effected by outlier.
   
clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':peakfit9demo ']);diary off
clear;peakfit9demo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':AsymmetricalAreaTest  ']);diary off
clear;AsymmetricalAreaTest

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':AsymmetricalAreaTest2  ']);diary off
clear;AsymmetricalAreaTest2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':FouFilterExampleScript  ']);diary off
clear;FouFilterExampleScript

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': MorseCode Demonstration']);diary off
clear;
MorseCode

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoSegmentedSharpenAreaRatio  ']);diary off
clear;DemoSegmentedSharpenAreaRatio

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoSegmentedSharpenAreaRatio2  ']);diary off
clear;DemoSegmentedSharpenAreaRatio2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SharpenedGaussianDemo ']);diary off
clear;SharpenedGaussianDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoSegmentedSharpen  ']);diary off
clear;DemoSegmentedSharpen

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoSegmentedSharpen2  ']);diary off
clear;DemoSegmentedSharpen2

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PeakCalibrationCurve  ']);diary off
clear;PeakCalibrationCurve

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoSegmentedSmooth  ']);diary off
DemoSegmentedSmooth

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DemoPeakfit  ']);diary off
clear;DemoPeakfit

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SharpenedGaussianDemo4terms  ']);diary off
clear;SharpenedGaussianDemo4terms

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SharpenedLorentzianDemo  ']);diary off
clear;SharpenedLorentzianDemo

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':SharpenedLorentzianDemo4terms  ']);diary off
clear;SharpenedLorentzianDemo4terms

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':GaussVsExpGauss  ']);diary off
clear;GaussVsExpGauss

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DeconvDemo4  ']);diary off
clear;DeconvDemo4

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DeconvDemo7  ']);diary off
clear;DeconvDemo7

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':DeconvDemo6  ']);diary off
clear;DeconvDemo6

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PowerMethodGaussian  ']);diary off
clear;PowerMethodGaussian

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PowerMethodLorentzian  ']);diary off
clear;PowerMethodLorentzian

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':PowerLawTest  ']);diary off
clear;PowerLawTest

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': MeasuringWidth ']);diary off
clear;MeasuringWidth

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': range ']);diary off
clear;range(1:99)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': ipf ']);diary off
clear;x=[0:.005:1];y=humps(x);ipf(x,y)

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': TFitCalCurveAbs ']);diary off
clear;TFitCalCurveAbs

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': makeodd ']);diary off
clear; makeodd([1.1 2 3 4.8 5 6 7.7 8 9])

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ': halfwidth ']);diary off
clear;
x=-5:.01:5;
y=sin(x).^2;
FWHM=halfwidth(x,y,0) % Area of central peak
FWHM=halfwidth(x,y,1.5) % Area of positive side peak

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':  ']);diary off
clear;

clear;load fn;fn=fn+1;save fn;diary on;disp([num2str(fn) ':  ']);diary off
clear;
clear;diary on;
endtime=cputime;
load TestStartTime
TotalElapsedTime=(endtime-TestStartTime)./60
diary off
delete fn.mat
disp('--------------------------------------------------------')