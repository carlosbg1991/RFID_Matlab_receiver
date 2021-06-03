increment=1;
x=[1:increment:4000];

% For each simulated peak, compute the amplitude, position, and width
amp=round(5.*randn(1,38));  % Amplitudes of the peaks  (Change if desired)
pos=[200:100:3900];   % Positions of the peaks (Change if desired)
wid=40.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=0.2; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.2  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
  end
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix