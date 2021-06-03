% Demonstration script for iSignal function. It generates a test signal
% consisting of several peaks, adds random noise, runs iSignal, then prints out 
% the actual peak positions, heights, widths and areas. 
% Each time you run it you get a different set of peaks. You can easily
% evaluate the accuracy of the measurements because the actual peak
% parameter values in this simulation are always integers.
%   T. C. O'Haver, July 2011
increment=1;
x=[1:increment:4000];

% For each simulated peak, compute the amplitude, position, and width
amp=round(5.*randn(1,38));  % Amplitudes of the peaks  (Change if desired)
pos=[200:100:3900];   % Positions of the peaks (Change if desired)
wid=60.*ones(size(pos));   % Widths of the peaks (Change if desired)

% Optionally add noise and background (Change if desired) 
Noise=.01; % Amount of random noise added to the signal. (Change if desired) 
background=10; % Strength of bnackground added to the signal. (Change if desired)

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.2,  % Keep only those peaks above a certain amplitude
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
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
isignal(demodata);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks




