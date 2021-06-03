% Here's an example of real-time smoothing that you can run on your
% computer in Matlab or Octave without any additional hardware.  It uses
% the mouse pointer to simulate an analog-to-digital converter. When you
% run this script, it will display a graph. To simulate a signal coming
% from an instrument or sensor, position your mouse pointer along the y
% (vertical) axis and click to enter 20 data points as you move the mouse
% pointer up and down the y axis simulating a signal. (If you wish you can
% change 'maxx' in line 17 and SmoothWidth in line 20).
% 
% To do this with real time data from a sensor, replace lines 25 and 26
% with the code that acquires one data point from your sensor.
% Tom O'Haver (toh@umd.edu) 2018
format compact
format short g
clf
clear
maxx=30; % Sets the maxumum number of data points
axis([0 maxx -10 10]);
hold on
SmoothWidth=7; % Smooth width (odd integer)
% SmoothWidth=makeodd(SmoothWidth); % Optional to insure oddness.
y=zeros(1,maxx); % y is vector of simulated data points
sy=zeros(1,maxx); % sy is vector of smoothed data points
for n=1:maxx
    [clickX,clickY] = ginput(1); % Wait for a mouse click
    y(n)=clickY; % simulate a data point by the y position of the click on the graph
    FirstPoint=n-SmoothWidth;
    if FirstPoint<1;FirstPoint=1;end
    LastPoint=n-1;
    if LastPoint<1;LastPoint=1;end
    SmoothedPoint=mean(y(FirstPoint:LastPoint));
    sy(n)=SmoothedPoint;
    if n==1
        plot(1,sy(1)) % Plot the first simulated data point
        xlabel(['X      SmoothWidth= ' num2str(SmoothWidth)]);ylabel('Y') % Label axes
        title('Real Time Smoothed Derivative Demo.  Black=original   Red=Smoothed Derivative')
    else
        LineStart=n-1;
        if LineStart<1;LineStart=1;end
        LineEnd=n;
        plot([n-1 n],[y(LineStart) y(LineEnd)],'k') % Draw raw signal as black line
         dy=deriv1(sy);
        plot([n-1 n],[dy(LineStart) dy(LineEnd)],'r') % Draw smoothed signal as red line
    end  
end
hold off

% REQUIRED FUNCTION. Place this function the the Matlab path.
% function d=deriv1(a)
% % First derivative of vector using simple adjacent differences.
% % Example: deriv([1 1 1 2 3 4]) returns [0 0 0 1 1 1]
% %  T. C. O'Haver, 2013.
% n=length(a);
% d=zeros(size(a));
% d(1)=a(2)-a(1);
% for j = 2:n;
%   d(j)=a(j)-a(j-1);
% end


