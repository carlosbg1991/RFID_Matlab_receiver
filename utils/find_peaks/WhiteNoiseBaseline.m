%  The script WhiteNoiseBaseline.m simulates a Gaussian peak with randomly
%  variable position and width, plus white noise with a
%  signal-to-noise ratio is 15, which is measured by least-squares curve
%  fitting methods using peakfit.m with different methods of baseline
%  correction:
% (a) a single-component model (shape 1) with autozero set to 1 (meaning a
% linear baseline is first interpolated from the edges of the data segment
% and subtracted from the signal): peakfit([x;y],0,0,1,1,0,10,1); 
% (b) a 2-component model, the first being a Gaussian (shape 1) and the
% second a linear slope (shape 26). with autozero set to 1:
% peakfit([x;y],0,0,2,[1 26],[0 0],10,0). The results are about the same
% for both methods, but the second method works better if the peak is near
% the edges of the data range.
% Compare to RandomWalkBaseline.m
PeakHeight=1;
PeakPosition=100+20.*randn();
PeakWidth=50+10.*randn();
x=1:200;
noise=.06.*randn(size(x));
baseline=noise; %  white noise
y=PeakHeight.*gaussian(x,PeakPosition,PeakWidth)+baseline;
SignalToNoiseRatio=PeakHeight/std(baseline)
figure(1)
% Single Gaussian peak fit, baseline mode 1
[FitResults1,GOF1]=peakfit([x;y],0,0,1,1,0,10,0,1,0,1,1);
subplot(2,1,1)
title('Gaussian on random walk: Single peak fit, baseline mode 1')
figure(2)
% Two-peak fit, one Gaussian + linear slope, baseline mode 0
[FitResults2,GOF2]=peakfit([x;y],0,0,2,[1 26],[0 0],10,0,0,0,1,1);
subplot(2,1,1)
title('Gaussian on random walk: Two-peak fit, one Gaussian + linear slope, baseline mode 0')
PositionError1=abs(100.*(PeakPosition-FitResults1(1,2))./PeakPosition);
HeightError1=abs(100.*(PeakHeight-FitResults1(1,3))./PeakHeight);
WidthError1=abs(100.*(PeakWidth-FitResults1(1,4))./PeakWidth);
PositionError2=abs(100.*(PeakPosition-FitResults2(val2ind(FitResults2(:,2),AvgPeakPosition),2))./PeakPosition);
HeightError2=abs(100.*(PeakHeight-FitResults2(val2ind(FitResults2(:,2),AvgPeakPosition),3))./PeakHeight);
WidthError2=abs(100.*(PeakWidth- FitResults2(val2ind(FitResults2(:,2),AvgPeakPosition),4))./PeakWidth);

disp('       Position Error  Height Error  Width Error')
disp(['Method 1:  ' num2str([PositionError1 HeightError1 WidthError1] )] )
disp(['Method 2:  ' num2str([PositionError2 HeightError2 WidthError2] )] )