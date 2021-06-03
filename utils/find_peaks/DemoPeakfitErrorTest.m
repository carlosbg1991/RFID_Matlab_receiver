% Script for demonstrating the peak parameters errors
% involved in using peakfit.m, with built-in signal generator.  
% See http://terpconnect.umd.edu/~toh/spectrum/CurveFittingC.html
% and http://terpconnect.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
% T. C. O'Haver (toh@umd.edu).  May 11, 2009
 
clear
close
format long g
format compact
warning off all

% Initial settings of users-adjustable parameters
n=300;  % maximum x-value
increment=1; % Difference between adjacent x values
x=[1:increment:n];  % Creates x-vector
amp=[1 2 3];  % Amplitudes of the peaks (one entry for each peak)
pos=[100 130 200];   % Positions of the peaks (one entry for each peak)
wid=[30 40 50];   % Widths of the peaks (one entry for each peak)
noise=.01; % Amplitude of random noise added to signal
NumRepeats=10; % Number of repeat samples with independent noise (minimun=2)
NumTrials=1; % Number of trail fits per sample
PeakShape=1; % Gaussian=1, Lorentzian=2, Logistic=3, Pearson=4
% "extra" determines the shape of the Pearson function. When extra is
% 1.0, the shape is Lorentzian. When extra is large (>20), the shape 
% approaches Gaussian. At small values (<1), the shape is a pointy "cusp".
extra=2; 
center=160; % x-value in center of fitted range
window=270; % range of x-values in fitted range
startvector=[98 28 129 41 190 49]; % First guess vector [pos1 wid1 pos2 wid2 ...]

% Generate simulated signal peaks
% A = matrix containing one of the unit-amplidude functions in each of its rows
clear A
A=zeros(length(pos),length(x));
for k=1:length(pos),
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
  end % switch
  area(k)=amp(k).*trapz(x,A(k,:)); % Computes peak areas via trapaziod method
end % for
TrueY=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
% TrueY=TrueY+1+(1./(x+10)); % Adds sloping curved baseline

% Arrange the peak parameters of the simulated signal into a matrix
ColumnLabels='           Peak#   Position     Height     Width       Area';
NumPeaks=length(pos);
for m=1:NumPeaks,
  if m==1,  
     ActualParameters=[m pos(m) amp(m) wid(m) area(m) ];, 
  else
     ActualParameters=[ActualParameters; [m pos(m) amp(m) wid(m) area(m) ]];
  end
end
   
% Generate and fit NumRepeats separate samples of the simulated signal with noise.
AllFitResults=zeros(NumPeaks,5,NumRepeats);
AllFittingErrors=zeros(1,NumRepeats);
AllPercentParameterErrors=zeros(NumPeaks,5,NumRepeats);
for repeat=1:NumRepeats,
   %  Add random noise to simulated signal
   y = TrueY + noise .* randn(size(x));
   % Execute the peakfit script
      startnow=startvector;
   delta=window/50;
   for k=1:2*NumPeaks,
        startnow(k)=startvector(k)+randn*delta;
    end
   [FitResults,LowestError]=peakfit([x' y'],center+randn*delta,window+randn*delta,NumPeaks,PeakShape,extra,NumTrials,startnow); 
   AllFitResults(:,:,repeat)=FitResults;
   AllFittingErrors(repeat)=LowestError(1);
   AllPercentParameterErrors(:,:,repeat)=100.*(FitResults-ActualParameters)./ActualParameters;
   % [AllFittingErrors(repeat) AllPercentParameterErrors(1,2:4,repeat)]
end

% Calculate the mean and relative standard deviaiton of parameters
for peak=1:NumPeaks,
      MeanParameters(peak,:)=mean(squeeze(AllFitResults(peak,:,:))');
      RSDparameters(peak,:)=100.*std(squeeze(AllFitResults(peak,:,:))')./mean(squeeze(AllFitResults(peak,:,:))');
end
  disp(['Noise amplitude: ' num2str(noise)])
  disp(['Number of repeats samples: ' num2str(NumRepeats)])
  disp(['Number of trial fits per sample: ' num2str(NumTrials)])
  disp('Peak#    Position        Height       Width       Area');
  disp('Actual parameters');disp(num2str(ActualParameters))
  disp('Average measured parameters');disp(num2str(MeanParameters))
  RSDparameters(:,1)=ActualParameters(:,1);
  disp('Percent RSD of measured parameters');disp(num2str(RSDparameters))
  disp(['Average percent fitting error: ' num2str(mean(AllFittingErrors))])
  PercentParameterErrors=100.*(MeanParameters-ActualParameters)./ActualParameters;
  PercentParameterErrors(:,1)=ActualParameters(:,1);
  disp('Percent error of measured parameters');disp(num2str(PercentParameterErrors))
  AverageParameterErrors=mean(abs(PercentParameterErrors(:,2:5)))
  disp(['Overall Average Parameter Error: ' num2str(mean(AverageParameterErrors)) ])
  xa=AllFittingErrors;
pe=squeeze(AllPercentParameterErrors(1,2,:));
he=squeeze(AllPercentParameterErrors(1,3,:));
we=squeeze(AllPercentParameterErrors(1,4,:));
fp=polyval(polyfit(xa,pe',1),xa);
fh=polyval(polyfit(xa,he',1),xa);
fw=polyval(polyfit(xa,we',1),xa);
figure(2);plot(xa,pe,'.r',xa,fp,'r',xa,he,'.b',xa,fh,'b',xa,we,'.g',xa,fw,'g')
xlabel('Percent Fitting Error')
ylabel('Percent Parameter Error')
title('Peak #1: RED = position   BLUE = height    GREEN = width')

pe=squeeze(AllPercentParameterErrors(2,2,:));
he=squeeze(AllPercentParameterErrors(2,3,:));
we=squeeze(AllPercentParameterErrors(2,4,:));
fp=polyval(polyfit(xa,pe',1),xa);
fh=polyval(polyfit(xa,he',1),xa);
fw=polyval(polyfit(xa,we',1),xa);
figure(3);plot(xa,pe,'.r',xa,fp,'r',xa,he,'.b',xa,fh,'b',xa,we,'.g',xa,fw,'g')
xlabel('Percent Fitting Error')
ylabel('Percent Parameter Error')
title('Peak #2: RED = position   BLUE = height    GREEN = width')

pe=squeeze(AllPercentParameterErrors(3,2,:));
he=squeeze(AllPercentParameterErrors(3,3,:));
we=squeeze(AllPercentParameterErrors(3,4,:));
fp=polyval(polyfit(xa,pe',1),xa);
fh=polyval(polyfit(xa,he',1),xa);
fw=polyval(polyfit(xa,we',1),xa);
figure(4);plot(xa,pe,'.r',xa,fp,'r',xa,he,'.b',xa,fh,'b',xa,we,'.g',xa,fw,'g')
xlabel('Percent Fitting Error')
ylabel('Percent Parameter Error')
title('Peak #3: RED = position   BLUE = height    GREEN = width')