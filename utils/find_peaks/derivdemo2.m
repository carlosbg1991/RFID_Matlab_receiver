% Matlab/Octave code example

function derivdemo2


% Create an independent variable ("x-axis")
increment=1;
x=[1:increment:300];

% Create a signal consisting of two Gaussian peaks, one twice as wide as the other:
% For each simulated peak, specify the amplitude, position, and width
amp=[1 1 ];  % Amplitudes of the peaks  (Change if desired)
pos=[100 200];   % Positions of the peaks (Change if desired)
wid=[20 40];   % Widths of the peaks (Change if desired)
Noise=.0; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix


% Plot the signal and its 1st, 2nd, and 3rd derivatives in the four quadrants of the plot window
subplot(2,2,1)
plot(x,y)
subplot(2,2,2)
plot(x,deriv(y))
subplot(2,2,3)
plot(x,deriv2(y))
subplot(2,2,4)
plot(x,deriv3(y))


% Internal sub-functions

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
%  Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
%  plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);

function d=deriv(a)
% First derivative of vector using 2-point central difference.
% Example: deriv([1 1 1 2 3 4]) yeilds [0 0 .5 1 1 1]
%  T. C. O'Haver, 1988.
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end

function d=deriv2(a)
% Second derivative of vector using 3-point central difference.
%  T. C. O'Haver, 2006.
n=length(a);
for j = 2:n-1;
  d(j)=a(j+1) - 2.*a(j) + a(j-1);
end
d(1)=d(2);
d(n)=d(n-1);

function d3=deriv3(a)
% Third derivative of vector a
%  T. C. O'Haver, 2008.
n=length(a);
for j = 3:n-2;
  d3(j)=a(j+2) - 2.*a(j+1) + 2.*a(j-1) - a(j-2);
end
d3(1:2)=d3(3);
d3(n-1:n)=d3(n-2);