% Demoautopeaksplot 
% A simple self-contained demonstration of the segmented peak finding
% function findpeaksSb function applied to noisy synthetic data set
% consisting of a random number of variable-width Gaussian peaks whose
% widths increase by a factor of 25-fold and that are superimposed in a
% curved baseline with random white noise that increases gradually; four
% segments are used, changing the peak detection and curve fitting values
% so that all the peaks are measured accurately. Each time you run this, a
% different set of peaks is generated. Calls the findpeaksSb and halfwidth.m
% functions, which must be in the Matlab path. See
% http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% Tom O'Haver (toh@umd.edu). Version 1 November 2016
% You can change the signal characteristics in lines 20-24
format short g
format compact
clear PercentError Measuredpeaks
increment=2;

% Initial values of variable parameters
SlopeThreshold=[0.001 .0002 .0001 .00005]; % 
AmpThreshold=[.00001 .001 .01 .1];
SmoothWidth=[5 30 70 80];  % SmoothWidth should be roughly equal the peak width (in points)
FitWidth=[10 12 15 40]; % FitWidth should be roughly equal to the peak widths (in points)
% Autozero sets the baseline correction mode in the last argument: autozero=0
% (default) does not subtract baseline from data segment;. autozero=1
% interpolates a linear baseline from the edges of the data segment and
% subtracts it from the signal (assumes that the peak returns to the
% baseline at the edges of the signal); autozero=2,  like mode 1 except
% that it computes a quadratic curved baseline; autozero=3 compensates for
% a flat baseline without reference to the signal itself (does not require
% that the signal return to the baseline at the edges of the signal, as
% does modes 1 and 2).
extra=0;
NumTrials=1;
autozero=3; 
x=1:increment:16000;

% For each simulated peak, compute the amplitude, position, and width
pos=500:400:16000;   % Positions of the peaks (Change if desired)
amp=round(2.*randn(1,length(pos))); % Amplitudes of the peaks  (Change if desired)
wid=30.*ones(size(pos)).*(pos/2000); % Increasing widths of the peaks (Change if desired)
Noise=.01; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.01,  % Keep only those peaks above a certain amplitude
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
y=z+Noise.*randn(size(z));  % Adds constant random noise
y=z+x/2000.*Noise.*randn(size(z));  % Adds constant random noise
% y=z+Noise.*sqrtnoise(z);  % Adds signal-dependent random noise
% demodata=[x' y']; % Assembles x and y vectors into data matrix

figure(1);plot(x,y,'r')  % Graph the signal in red
title('autopeaksplot demo. Detected peaks are numbered. Peak table is printed in Command Window')
% ActualPeaks

% Label the x-axis with the parameter values

% Find and measure the peaks
tic;
Measuredpeaks=autopeaks(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
ElapsedTime=toc;
PeaksPerSecond=length(Measuredpeaks)/ElapsedTime;
% Measure the halfwidth of each peak
NumPeaks=max(Measuredpeaks(:,1));
for n=1:NumPeaks;
    Measuredpeaks(n,4)=halfwidth(x,y,Measuredpeaks(n,2));
end

% Display results
disp('-------------------autopeaksplot.m demo--------------------')
disp(['SlopeThreshold = ' num2str(SlopeThreshold) ] )
disp(['AmpThreshold = ' num2str(AmpThreshold) ] )
disp(['SmoothWidth = ' num2str(SmoothWidth) ] )
disp(['FitWidth = ' num2str(FitWidth) ] )
disp(['Speed = ' num2str(round(PeaksPerSecond)) ' Peaks Per Second' ] )
% disp('         Peak #     Position      Height      Width       Area      % Fitting error   R2')
% MeasuredPeaks=Measuredpeaks(:,1:7)  % Display table of peaks
% calculate top of peaks including background for peak number labeling
clear peaky
SizeMeasuredpeaks=size(Measuredpeaks);NumPeaksFound=SizeMeasuredpeaks(1);
for peak=1:NumPeaksFound,peaky(peak)=y(val2ind(x,Measuredpeaks(peak,2)));end
% Number the peaks found on the graph
figure(1);text(Measuredpeaks(:, 2),peaky',num2str(Measuredpeaks(:,1)))
% xlabel(['SlopeThresh.=[' num2str(SlopeThreshold) ']    AmpThresh.=[' num2str(AmpThreshold) ']    SmoothWidth=[' num2str(SmoothWidth) ']    FitWidth=[' num2str(FitWidth) '] ' ])
for peak=1:NumPeaksFound,
    for param=1:6
        TrueValue=ActualPeaks(val2ind(ActualPeaks(:,param),Measuredpeaks(peak,param)),param);
        Error=100*(TrueValue-Measuredpeaks(peak,param))./TrueValue;
        PercentError(peak,param)=Error(1);
    end
end
%disp('Relative Percent Errors')
% disp('     Position       Height        Width        Area')
% disp(PercentError(:,2:5))
 disp('Average Relative Percent Errors')
  disp('     Position       Height        Width        Area')
 disp(sqrt(mean(PercentError(:,2:5).^2)))
% figure(2);ipf(x,y);
