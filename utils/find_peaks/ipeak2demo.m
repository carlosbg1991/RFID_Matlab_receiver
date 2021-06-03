% Demonstration script for iPeak2 function.  This is a demo script for 
% ipeak that demonstrates the operation and the typical measurement 
% errors of ipeak.  It generates a test signal consisting of several 
% Gaussian peaks, adds random noise, runs ipeak, then prints out the 
% the actual and the measured peak positions, heights, widths and areas. 
% Each time you run it you get a different set of peaks. You can easily
% evaluate the accuracy of the measurements because the actual peak
% parameter values in this simulation are always integers.
%   T. C. O'Haver, June 2011
increment=2;
x=[1:increment:4000];

% For each simulated peak, compute the amplitude, position, and width
amp=round(5.*randn(1,38));  % Amplitudes of the peaks  (Change if desired)
pos=[200:100:3900];   % Positions of the peaks (Change if desired)
wid=70.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=.05; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.2,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x positions
      % A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))).^2); % Gaussian peaks
      A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles peak parameters into ActualPeaks matrix
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) amp(k).*wid(k)]; 
      p=p+1;
  end; 
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=ipeak2(demodata,0,0.5,0.0015,25,20);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks


