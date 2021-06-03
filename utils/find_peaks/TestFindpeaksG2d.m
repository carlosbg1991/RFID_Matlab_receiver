% Demonstration of findpeaks2d.m, which shows that this function can locate
% peaks resulting in 'shoulders' that do not produce a distinct maximum in
% the original signal. Detects peaks in the negative of the smoothed second
% derivative of the signal (shown as the dotted line in the figure).
% Requires gaussian.m, findpeaksG.m, findpeaksG2d.m, fastsmooth.m, and
% peakfit.m in the path. Download from http://tinyurl.com/cey8rwh
format compact
figure(1);
x=1:.1:50;
noise=.1; % << change: Try larger or smaller amounts of noise
AmpThreshold=0.4;
SlopeThreshold=1e-7;
smoothwidth=25;
y1=gaussian(x,20,10); % Peak 1: Position=20 Height=1 Width=10
y2=.5.*gaussian(x,30,10); % Peak 2: Position=30 Height=0.5 Width=10
n=noise.*randn(size(x)); % random white noise
y=y1+y2+n; % Simulated signal contains two peaks plus noise
clf;
sy=fastsmooth(y,smoothwidth,3); % Pre-smoothing before calling findpeaks
ssd=-2000.*fastsmooth(deriv2(sy),smoothwidth,3); % -Smoothed second derivative
plot(x,y1,x,y2,x,y,x,ssd,'-.');
disp(' ')
noise=noise
disp('    Peak     Position   Height     Width     Area')
PfindpeaksG=findpeaksG(x,sy,SlopeThreshold,AmpThreshold,smoothwidth,25,3)

PfindpeaksG2d=findpeaksG2d(x,sy,SlopeThreshold,AmpThreshold,smoothwidth,25,3)
xlabel('Blue= Peak 1    Green= Peak 2    Red= Peak 1 + Peak 2    Dotted= -smoothed second derivative');
title('Comparison of findpeaksG.m and findpeaksG2d.m for locating a weak side shoulder')
disp('findpeaksG detects only one peak, but findpeaksG2d detects both. ')
disp(' ')
disp('If the findpeaksG2d results are not accurate enough, you can use its ')
disp('results as the "start" value for peakfit.m, which takes longer but gives')
disp('more accurate results, especially for width and area:')
% Comparison to peakfit
sizeP=size(PfindpeaksG2d);
NumPeaks=sizeP(1);
start=[PfindpeaksG2d(1,2) PfindpeaksG2d(1,4) PfindpeaksG2d(2,2) PfindpeaksG2d(2,4)];
figure(2)
disp('PeakfitResults=peakfit([x;y],0,0,NumPeaks,0,10,start)')
PeakfitResults=peakfit([x;y],0,0,NumPeaks,0,10,start)