function peakfitdemob
% Demonstration function for peakfit function. It generates a test signal
% consisting of 3 peaks on a curved baseline, then runs peakfit function,
% treating the baseline as an 4th peak whose peak position is negative.
%   T. C. O'Haver, December 2013
increment=.1;
x=[1:increment:500];

% For each simulated peak, compute the amplitude, position, and width
% Peak 1 is the baseline
amp=[1 2 3 100];  % Amplitudes of the peaks  (Change if desired)
pos=[100 250 400 -200];   % Positions of the peaks (Change if desired)
wid=[50 50 50 700];   % Widths of the peaks (Change if desired)
Noise=.01; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplitude peak in each of its rows
A = zeros(length(pos),length(x));
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      A(k,:)=gaussian(x,pos(k),wid(k)); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
demodata=[x' y']'; % Assembles x and y vectors into data matrix

disp(' ')
disp('Measurement of three Gaussian peaks with heights 1, 2, and 3,')
disp('superimposed on a very strong broad Gaussian baseline.')
disp('The signal is fit with the peakfit function for four peaks, treating')
disp('the baseline as an additional peak whose peak position is negative.')
disp(' ')
disp(['Signal to noise ratio of smallest peak=' num2str(min(amp./Noise))])
% Call peakfit to fit 4 peaks, with the 4th one intended to fit the
% baseline. This requires a specified start value, because the default
% start assumes that all peaks are within the data range
%
clf
NumPeaks=length(amp);
PeakShape=1; % 1-Gaussian, 2=Lorentzian, etc.
start=[100 50 250 50 400 50 -200 700];
NumTrials=1; % Try something higher if 1 is unreliable
BaselineMode=0;
tic
PeakfitResults=peakfit(demodata,250,500,NumPeaks,PeakShape,0,NumTrials,start,BaselineMode);
toc
disp('      Peak 1        Peak 2       Peak 3       Baseline ')
disp(PeakfitResults(:,3)')
disp(' ')
disp('You can test the reliability of this method by changing')
disp('the peak parameters in lines 11, 12, and 13 and see if the')
disp('peakfit function will track the changes WITHOUT having to')
disp('change the start vector in liine 42.')
%
% Compare to the Classical Least Squares (CLS) technique assuminig that the
% peak positions and widths of all peaks, including the baseline are
% known beforehand.
%
disp(' ')
disp('Compare to the Classical Least Squares (CLS) technique:')
tic
CLSresults=cls(demodata(1,:),demodata(2,:),NumPeaks,PeakShape,pos,wid,0);
toc
CLSresults
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);

% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);