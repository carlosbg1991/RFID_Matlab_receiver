% Here's an example of real-time smoothing that you can run on your
% computer in Matlab or Octave without any additional hardware.  Lines 9-25
% generate a simulated signal and a smoothing vector. Lines 27-37 perform
% the actual point-by-point smoothing operation. (If you wish you can
% change 'maxx' in line 16 to control the total number of data points and
% SmoothWidth in line 17). Uncomment either line 24 (rectangular), 25
%  (triangular), or 26 (Gaussian), or write a smooth vector of your choice.
% To do this with real time data from a sensor, comment out line 25 and
% replace line 28 with the code that acquires one data point from your
% sensor. This version acquires, smooths, and stores the smoothed data in
% sy in real time, then plots the data only after data acquisition.
% Needs triangle and rectanglepulse functions. Tom O'Haver (toh@umd.edu) 2018
format compact
format short g
clear
maxx=500; % Sets the maxumum number of data points
SmoothWidth=41; % Smooth width (odd integer)
% SmoothWidth=makeodd(SmoothWidth); % Optional to insure oddness.
% Uncomment one of the following line to define desired smooth shape.
% SmoothVector=ones(1,SmoothWidth);SmoothVector=SmoothVector./(sum(SmoothVector));
SmoothVector=triangle(1:SmoothWidth,SmoothWidth/2,SmoothWidth/2);SmoothVector=SmoothVector./(sum(SmoothVector));
% SmoothVector=gaussian(1:SmoothWidth,SmoothWidth/2,SmoothWidth/2);SmoothVector=SmoothVector./(sum(SmoothVector));
y=zeros(1,maxx); % y is vector of simulated data points
sy=zeros(1,maxx); % sy is vector of smoothed data points
RawData=rectanglepulse(1:maxx,maxx/2,maxx/4)+.02.*randn(size(y)); % Generate sumulated data
tic
for n=1:maxx
    y(n)=RawData(n);  % <<< replace this line with call to your data acquisition harware
    FirstPoint=n-SmoothWidth;
    if FirstPoint<1;FirstPoint=1;end
    LastPoint=n-1;
    if LastPoint<1;LastPoint=1;end
    if n>SmoothWidth
        SmoothedPoint=sum(SmoothVector.*(y(FirstPoint:LastPoint)));
        sy(n)=SmoothedPoint;
    end
end
ElapsedTime=toc
TimePerPoint=ElapsedTime./maxx
plot(1:maxx,sy) % Plot the entire smoothed signal
xlabel(['X      SmoothWidth= ' num2str(SmoothWidth)]);ylabel('Y') % Label axes
title('Real Time Smooth Demo.  Black=original   Red=Smoothed')

% On a standard desktop PC (Intel Core i5 3 Ghz) running Windows 10 home,
% the elapsed time per data point is 2 microseconds (without plotting)
% and 1.4 milliseconds per point with point-by-point plotting (lines 40-50).


