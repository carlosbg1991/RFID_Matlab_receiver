% Here's an example of real-time smoothed peak detection that you can run
% on your computer in Matlab or Octave without any additional hardware.  It
% uses the mouse pointer to simulate an analog-to-digital converter. When
% you run this script, it will display a graph. To simulate a signal coming
% from an instrument or sensor, position your mouse pointer along the y
% (vertical) axis and click to enter 20 data points as you move the mouse
% pointer up and down the y axis simulating a signal. (If you wish you can
% change 'maxn' in line 17 and SmoothWidth in line 20).
%
% To do this with real time data from a sensor, replace lines 25 and 26
% with the code that acquires one data point from your sensor.
% Requires triangle and gaussian functions. Tom O'Haver (toh@umd.edu) 2018
format compact
format short g
clf
clear
load DataMatrix2
maxn=3000; % Sets the maxumum number of data points
axis([0 maxn 0 10]);
hold on
SmoothWidth=5; % Smooth width (odd integer)
% SmoothWidth=makeodd(SmoothWidth); % Optional to insure oddness.
x=zeros(1,maxn); % x is vector of simulated data points x=DataMatrix2(:,1);
y=zeros(1,maxn); % y is vector of simulated data points y=DataMatrix2(:,2);
sy=zeros(1,maxn); % sy is vector of smoothed data points
for n=1:maxn
    x(n)=DataMatrix2(n,1);
    y(n)=DataMatrix2(n,2); % simulate a single data point from the data source
    FirstPoint=n-SmoothWidth;
    if FirstPoint<1;FirstPoint=1;end
    LastPoint=n-1;
    if LastPoint<1;LastPoint=1;end
    SmoothedPoint=mean(y(FirstPoint:LastPoint));
    sy(n)=SmoothedPoint;
    if n==1
        plot(1,sy(1)) % Plot the first simulated data point
        xlabel(['X      SmoothWidth= ' num2str(SmoothWidth)]);ylabel('Y') % Label axes
        title('Real Time Smoothed Peak Detection Demo.  Black=original   Red=Smoothed')
    else
        LineStart=n-1;
        if LineStart<1;LineStart=1;end
        LineEnd=n;
        plot([n-1 n],[sy(LineStart) sy(LineEnd)],'r') % Draw smoothed signal as red line
        plot([n-1 n],[y(LineStart) y(LineEnd)],'k') % Draw raw signal as black line
        if n>2 % Start peak detection when 3 points have been acquired.
            if sy(n-1)>sy(n-2) % If a point is greater that the the previous one
                if sy(n-1)>sy(n) % AND greater than the following one, register a peak.
                    disp(['Peak detected at x=' num2str(n-1) ' and y=' num2str(y(n-1))])
                    text(n-1,sy(n-1),' peak') % Label the peak on the graph
                end
            end
        end
    end
end
hold off


