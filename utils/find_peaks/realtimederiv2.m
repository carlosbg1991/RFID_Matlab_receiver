% Here's an example of real-time 2nd differentiation that you can run on
% your computer in Matlab or Octave without any additional hardware.  It
% uses the mouse pointer to simulate an analog-to-digital converter. When
% you run this script, it will display a graph. To simulate a signal coming
% from an instrument or sensor, position your mouse pointer along the y
% (vertical) axis and click to enter 20 data points as you move the mouse
% pointer up and down the y axis simulating one or more peaks. The original
% (raw) signal is plotted as a black line and its second derivative is
% plotted as a red line, (If you wish you can change 'maxx', the maximum
% number of data points, in line 16).
% 
% To do this with real time data from a sensor, replace lines 22 and 23
% with the code that acquires one data point from your sensor. Tom O'Haver
% (toh@umd.edu) 2018
clf
maxx=20;
axis([0 20 -10 10]);
hold on;
y=zeros(1,maxx);
dy=zeros(1,maxx);
for n=1:maxx-1
     [clickX,clickY] = ginput(1);
     y(n)=clickY;
    if n<3
        plot(1,y(1))
        xlabel('x');ylabel('y') % Label axes
        title('Real Time Second Derivative Demo')
    else
        plot([n-1 n],[y(n-1) y(n)],'k')
        dy(n)=-2.*y(n-1)+y(n-2)+y(n);
        plot([n-1 n],[dy(n-1) dy(n)],'r')
    end  
end
hold off

