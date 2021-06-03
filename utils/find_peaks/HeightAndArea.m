% HeightAndArea.m is a demonstration script that uses measurepeaks.m to
% measure the peaks in computer-generated signals consisting of a series of
% Gaussian peaks with gradually increasing widths that are superimposed on
% a curved baseline plus random white noise. It plots the signal (in Figure
% 1) and the individual peaks (in Figure 2) and compares the actual peak
% position, heights, and areas of each peak to those measured by
% measurepeaks.m using the absolute peak height, peak-valley difference,
% perpendicular drop, and tangent skim methods. Prints out a table of the
% relative percent difference between the actual and measured values for
% each peak and the average error for all peaks. To run, type "HeightAndArea".
% Note: measurepeaks.m function uses SmoothWidth only for peak detection;
% it performs measurements on the raw unsmoothed y data. In some cases it
% may be beneficial to smooth the y data yourself before calling
% measurepeaks.m, using any smooth function of your choice. (Uncomment 
% lines 94 to 122 to activate presmooting and compare to no smoothing).
% http://terpconnect.umd.edu/~toh/spectrum/Integration.html
% Tom O'Haver (toh@umd.edu).  Version 1.1 December 2016
% Requires measurepeaks.m, gaussian.m, val2ind.m in path.
clear
format short g
format compact
% Initial values of variable parameters
increment=5; % Interval between x values (smaller=more data points)
Noise=.05; % Amount of random noise added to the signal. (Change if desired) 
background=0; % Amplitude of broad background/baseline under the peaks
SlopeThreshold=.001; %  0.7*WidthPoints^-2, where WidthPoints is the number of data points in the half-width of the peak.  
AmpThreshold=.5; % Set just lower than the lowest peak height of interest
SmoothWidth=7;  % Normally roughly equal the peak width (in points)
FitWidth=11; % Normally roughly equal to the peak widths (in points)
%
x=1:increment:8000;
clf
% For each simulated peak, compute the amplitude, position, and width
pos=200:250:7500;   % Positions of the peaks (Change if desired)
% amp=ones(size(pos));
amp=round(5.*randn(1,length(pos))); % Amplitudes of the peaks  (Change if desired)
% wid=100.*ones(size(pos)).*sqrt(pos/20000); % Increasing widths of the peaks (Change if desired)
% wid=100.*ones(size(pos))./1.0646;
wid=300.*ones(size(pos)).*sqrt(pos/20000); % Increasing widths of the peaks (Change if desired)
% A = matrix containing one of the unit-amplitude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.9,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k)*wid(k) 0]; 
      p=p+1;
  end; 
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
% z=ones(size(pos))*A; % Substitute this to make all peak heights = 1.000
y=z+Noise.*randn(size(z));  % Adds constant random noise, stdev = Noise
% y=z+x/1000.*Noise.*randn(size(z));  % Adds constant random noise
% y=z+Noise.*sqrtnoise(z);  % Adds signal-dependent random noise
y=y+background*gaussian(x,-1000,8000); % Adds a broad background whose amplitude=background
% demodata=[x' y']; % Assembles x and y vectors into data matrix

Measured=measurepeaks(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,1);
sizeP=size(Measured);
NumPeaks=sizeP(1);


for peak=1:NumPeaks,
    PeakPosError(peak)=ActualPeaks(peak,2)-Measured(peak,2);
    PeakHeightError(peak)=ActualPeaks(peak,3)-Measured(peak,3);
    PeakPVError(peak)=ActualPeaks(peak,3)-Measured(peak,4);
    PeakPDAreaError(peak)=ActualPeaks(peak,5)-Measured(peak,5);
    PeakTSAreaError(peak)=ActualPeaks(peak,5)-Measured(peak,6);
end
 
 for peak=1:NumPeaks,
    PerecentPosError(peak)=-100.*(ActualPeaks(peak,2)-Measured(peak,2))./ActualPeaks(peak,2);
    PerecentHeightError(peak)=-100.*(ActualPeaks(peak,3)-Measured(peak,3))./ActualPeaks(peak,3);
    PerecentPVError(peak)=-100.*(ActualPeaks(peak,3)-Measured(peak,4))./ActualPeaks(peak,3);
    PerecentPDAreaError(peak)=-100.*(ActualPeaks(peak,5)-Measured(peak,5))./ActualPeaks(peak,5);
    PerecentTSAreaError(peak)=-100.*(ActualPeaks(peak,5)-Measured(peak,6))./ActualPeaks(peak,5);
    PercentErrors(peak,:)=[peak PerecentPosError(peak) PerecentHeightError(peak) PerecentPVError(peak) PerecentPDAreaError(peak) PerecentTSAreaError(peak)];
 end
 disp(' ')
 disp('Height and Area Relative Percent errors (Percent difference between actual and measured values)')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
 disp(PercentErrors)
disp(' ')
 disp('Average absolute relative percent errors of all peaks')
 disp('     Position     PeakMax      Peak-valley    Perp drop   Tan skim')
 disp([mean(abs(PerecentPosError)) mean(abs(PerecentHeightError))   mean(abs(PerecentPVError))  mean(abs(PerecentPDAreaError)) mean(abs(PerecentTSAreaError)) ] )
 
 
% % Pre-smooth data the data (uncomment the following lines to activate)
% PreSmoothWidth=3;
% y=SegmentedSmooth(y,PreSmoothWidth,3,1); 
% clear Measured
% Measured=measurepeaks(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,1);
% sizeP=size(Measured);
% NumPeaks=sizeP(1);
% 
% for peak=1:NumPeaks,
%     PeakPosError(peak)=ActualPeaks(peak,2)-Measured(peak,2);
%     PeakHeightError(peak)=ActualPeaks(peak,3)-Measured(peak,3);
%     PeakPVError(peak)=ActualPeaks(peak,3)-Measured(peak,4);
%     PeakPDAreaError(peak)=ActualPeaks(peak,5)-Measured(peak,5);
%     PeakTSAreaError(peak)=ActualPeaks(peak,5)-Measured(peak,6);
% end
%  
%  for peak=1:NumPeaks,
%     PerecentPosError(peak)=-100.*(ActualPeaks(peak,2)-Measured(peak,2))./ActualPeaks(peak,2);
%     PerecentHeightError(peak)=-100.*(ActualPeaks(peak,3)-Measured(peak,3))./ActualPeaks(peak,3);
%     PerecentPVError(peak)=-100.*(ActualPeaks(peak,3)-Measured(peak,4))./ActualPeaks(peak,3);
%     PerecentPDAreaError(peak)=-100.*(ActualPeaks(peak,5)-Measured(peak,5))./ActualPeaks(peak,5);
%     PerecentTSAreaError(peak)=-100.*(ActualPeaks(peak,5)-Measured(peak,6))./ActualPeaks(peak,5);
%     PercentErrors(peak,:)=[peak PerecentPosError(peak) PerecentHeightError(peak) PerecentPVError(peak) PerecentPDAreaError(peak) PerecentTSAreaError(peak)];
%  end
% 
%  disp('After pre-smoothing the y data')
%  disp('Average absolute relative percent errors of all peaks')
%  disp('     Position     PeakMax      Peak-valley    Perp drop   Tan skim')
%  disp([mean(abs(PerecentPosError)) mean(abs(PerecentHeightError))   mean(abs(PerecentPVError))  mean(abs(PerecentPDAreaError)) mean(abs(PerecentTSAreaError)) ] )
%  