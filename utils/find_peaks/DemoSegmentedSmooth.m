format short g
format compact
increment=10;
x=1:increment:8000;
clf
% For each simulated peak, compute the amplitude, position, and width
pos=300:400:7500;   % Positions of the peaks (Change if desired)
amp=round(5.*randn(1,length(pos))); % Amplitudes of the peaks  (Change if desired)
wid=100.*ones(size(pos)).*sqrt(pos/500); % Increasing widths of the peaks (Change if desired)
Noise=.5; % Amount of random noise added to the signal. (Change if desired) 
Baseline=2;
smoothwidths=[31 41 51];
SmoothType=3;
ends=1;
% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.1  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k)*wid(k)]; 
      p=p+1;
  end; 
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Adds constant random noise
% y=z+Noise.*sqrtnoise(z);  % Adds signal-dependent random noise
y=y+Baseline*gaussian(x,0,8000); % Optionally adds a broad background signal
% demodata=[x' y']; % Assembles x and y vectors into data matrix
sy=SegmentedSmooth(y,smoothwidths,SmoothType,ends);
NumSegments=length(smoothwidths);
figure(1);
subplot(1,2,1)
plot(x,y,'b')  % Graph the original signal in blue
thisaxis=axis;
xlabel('Original unsmoothed data')
title('SegmentedSmooth demo')
subplot(1,2,2)

plot(x,sy,'r')  % Graph the smoothed signal in red
axis(thisaxis);
xlabel('Smoothed by SegmentedSmooth.m')
title(['Number of segments = ' num2str(NumSegments) '    Smooth widths =  ' num2str(smoothwidths) ])
% disp('        Peak #     Position       Height        Width        Area')
% disp(ActualPeaks)