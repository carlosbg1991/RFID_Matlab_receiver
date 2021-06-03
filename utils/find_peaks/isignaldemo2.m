% Demonstration script for iSignal function. It generates a test signal
% consisting of 4 peaks, with heights 4, 3, 2 and 1, with equal widths,
% superimposed on a very strong curved baseline, plus added random white
% noise. The objective is to extract a measure that is proportional to the
% peak height but independent of the baseline strength. Suggested
% approaches: (a) Use automatic or manual baseline subtraction to remove
% the baseline; or (b) use differentiation (with smoothing) to suppress the
% baseline; or (c) use curve fitting (Shift-F) to measure peak height. 
% Each time you run it you get a different sample of random noise. 
%   T. C. O'Haver, March 2016

x=1:4000;

% For each simulated peak, compute the amplitude, position, and width
amp=[4 3 2 1];  % Amplitudes of the peaks  (Change if desired)
pos=[750 1500 2200 3000];   % Positions of the peaks (Change if desired)
wid=200.*ones(size(pos));   % Widths of the peaks (Change if desired)

% Optionally add noise and background (Change if desired) 
Noise=.1; % Amount of random noise added to the signal. (Change if desired) 
background=200; % Strength of background added to the signal. (Change if desired)

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.01,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.065.*amp(k).*wid(k)]; 
      p=p+1;
  end; 
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
y=y+background.*exp(-((x+2000)./(2402)).^2); % Optionally adds a broad background signal
DataMatrix=[x' y']; % Assembles x and y vectors into data matrix
figure(1)
isignal(DataMatrix);
disp('Demonstration script for iSignal function. It generates a test signal')
disp('consisting of 4 peaks, with heights 4, 3, 2 and 1, with equal widths,')
disp('superimposed on a very strong curved baseline, plus added random white')
disp('noise. The objective is to extract a measure that is proportional to the')
disp('peak height but independent of the baseline strength. Suggested')
disp('approaches: (a) Use automatic or manual baseline subtraction to remove')
disp('the baseline; or (b) use differentiation (with smoothing) to suppress the')
disp('baseline; or (c) use curve fitting (Shift-F) to measure peak height.')
disp('Each time you run this script, you get a different sample of random noise.')
disp(' ')

input('Press Enter to perform an automatic 3rd derivative calibration.')
SmoothMode=3; % Try other smooth widths also
SmoothWidth=150; % Try other smooth widths also
DerivativeMode=3; % Try other derivative modes also
P=isignal(DataMatrix,3000,500,SmoothMode,SmoothWidth,0,DerivativeMode,0,10,1000,0,0,0);
% Extract signal range from each peak of third derivative.
Xrange=500:1000; % x-axis range covering the maximum and minimum if the First peak
S1=max(P(Xrange))-min(P(Xrange)); % Peak-to-peak signal range over this interval
Xrange=1300:1700;% x-axis range covering the maximum and minimum if the Second peak
S2=max(P(Xrange))-min(P(Xrange));
Xrange=1900:2400;% x-axis range covering the maximum and minimum if the Third peak
S3=max(P(Xrange))-min(P(Xrange));
Xrange=2700:3300;% x-axis range covering the maximum and minimum if the Fourth peak
S4=max(P(Xrange))-min(P(Xrange));
figure(2)
clf
plotit(amp,[S1 S2 S3 S4],1,'ro');
xlabel('Actual peak heights')
ylabel('Peak-to-peak amplitude of derivatives')
title('Calibration curve based on the derivative of the signal')
