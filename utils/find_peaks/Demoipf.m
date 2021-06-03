% Demonstration script for ipf.m, with built-in signal generator.  
% You can change the character of the signal in lines 11-117.
% T. C. O'Haver (toh@umd.edu).  Version 3: August 23, 2008. 

close
format short g
format compact
warning off all

% Generate simulated signal
n=2000;  % maximum x-value
increment=2; % Difference between adjacent x values
x=[1:increment:n];
amp=[1 2 3 1 4 1 3 2 1 2 1];  % Amplitudes of the peaks (one entry for each peak)
pos=[100 130 200 400 600 630 800 850 900 1200 1700];   % Positions of the peaks (one entry for each peak)
wid=[30 40 50 30 20 20 30 40 50 100 200];   % Widths of the peaks (one entry for each peak)
noise=.1;
% A = matrix containing one of the unit-amplidude functions in each of its rows
clear A
A=zeros(length(pos),length(x));
for k=1:length(pos),
  A(k,:)=exp(-((x-pos(k))./(0.6005612.*wid(k))) .^2);  % Gaussian function
  area(k)=amp(k).*trapz(x,A(k,:)); % Computes peak areas via trapaziod method
end
TrueY=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
% TrueY=TrueY+lorentzian(x,0,1000); % Adds sloping curved baseline 
y = TrueY + noise .* randn(size(x));

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

% Display the peak parameters of the simulated signal 
ColumnLabels
ActualParameters

% Execute the Interactive Peak Fitter script
ipf(x,y);