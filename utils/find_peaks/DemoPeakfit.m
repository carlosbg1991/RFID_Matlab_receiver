% Demonstration script (for peakfit.m, version 3.1). Generates an
% overlapping peak signal, adds noise, fits it with peakfit.m, then
% computes the accuracy and precision of peak paramater measurements. 
% You can change any of the initial values in lines 11-24. T. C. O'Haver 
% (toh@umd.edu) Ver. 7: November, 2017.
format short g
format compact
warning off all
disp('Computing......... may take several seconds.')
% >>> Initial settings of user-adjustable parameters <<<
n=10;  % maximum x-value
increment=.1; % Difference between adjacent x values
x=1:increment:n;  % Creates x-vector
amp=[1 .5];  % Amplitude  s of the peaks (one entry for each peak)
pos=[4.8 6.2];   % Positions of the peaks (in increasing order)
wid=[1.666 1.666];   % Widths of the peaks (full width at half-maximum)
noise=.02; % Standard deviation of random noise added to signal
NumRepeats=10; % Number of repeat samples with independent noise
NumTrials=3; % Number of trail fits per sample
PeakShape=1; % Gaussian=1, Lorentzian=2, etc.
extra=2; 
center=5; % x-value in center of fitted range
window=10; % range of x-values in fitted range
startvector=[4 1 7 1]; % First guess vector [pos1 wid1 pos2 wid2 ...] Must be two values for each peak
baselinemode=0; % 0=no baseline correction; 1=linear baseline correction

% Generate simulated signal peaks
% A = matrix containing one of the unit-amplidude functions in each of its rows
clear A
A=zeros(length(pos),length(x));
for k=1:length(pos)
  switch PeakShape
  case 1
    A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))) .^2);  % Gaussian function
  case 2
    A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2); % Lorentzian function
  case 3
    L = exp(-((x-pos(k))/(.477.*wid(k))) .^2);
    A(k,:)= (2.*L)./(1+L);
  case 4
      % pearson function (shape determined by value of extra in line 21
      A(k,:)=ones(size(x))./(1+((x-pos(k))./((0.5.^(2/extra)).*4.62.*wid(k))).^2).^extra;
  case 6
    A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))) .^2);  % Gaussian function
  case 7
    A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Laorentzian function
  end % switch
  area(k)=amp(k).*trapz(x,A(k,:)); % Computes peak areas via trapaziod method
end % for
TrueY=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
% To test the autozero baseline subtraction, un-comment next line.
% TrueY=TrueY+ones(size(x))./(1+((x+100)./(300)).^2); % Adds sloping curved baseline

% Generate and fit NumRepeats separate samples of the simulated signal with noise.
NumPeaks=length(amp);
AllFitResults=zeros(NumPeaks,5,NumRepeats);
AllFittingErrors=zeros(1,NumRepeats);
AllPercentParameterErrors=zeros(NumPeaks,5,NumRepeats);
for repeat=1:NumRepeats
   %  Add different white noise sample to each repeat of simulated signal
   y = TrueY + noise .* randn(size(x));
   % Execute the peakfit function
   [FitResults,LowestError]=peakfit([x' y'],center,window,NumPeaks,PeakShape,extra,NumTrials,startvector,baselinemode); 
   FitResults=sortrows(FitResults,2); % Sort results by increasing peak positions (column 2)
   Heights1(repeat)=FitResults(1,3);
   Positions1(repeat)=FitResults(1,2);
   Widths1(repeat)=FitResults(1,4);
   Areas1(repeat)=FitResults(1,5);
   Heights2(repeat)=FitResults(2,3);
   Positions2(repeat)=FitResults(2,2);
   Widths2(repeat)=FitResults(2,4);
   Areas2(repeat)=FitResults(2,5);
end

% Calculate the mean and relative standard deviaiton of peak parameters
% measured by curve fitting. 
disp('Mean Fit Results')
disp('            Position     Height      Width       Area')
disp(['Peak 1       '  num2str([mean(Positions1) mean(Heights1) mean(Widths1) mean(Areas1)])  ] )
disp(['Peak 2       '  num2str([mean(Positions2) mean(Heights2) mean(Widths2) mean(Areas2)])  ] )
disp(' ') 
disp('Percent relative standard deviations')
disp('            Position     Height      Width       Area')
disp(['Peak 1       '  num2str([rsd(Positions1).*100 rsd(Heights1).*100 rsd(Widths1).*100 rsd(Areas1).*100])  ] )
disp(['Peak 2       '  num2str([rsd(Positions2).*100  rsd(Heights2).*100 rsd(Widths2).*100 rsd(Areas2).*100])  ] )

% Required functions: peakfit.m, rsd.m (shown below)

% function relstddev=rsd(x)
% % Relative standard deviation of vector x
% relstddev=std(x)./mean(x);
