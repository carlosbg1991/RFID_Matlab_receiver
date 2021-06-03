% Here's an example of real-time signal plotting that you can run on your
% computer in Matlab or Octave without any additional hardware.  It uses
% the mouse pointer to simulate an analog-to-digital converter. When you
% run this script, it will display a graph. To simulate a signal coming
% from an instrument or sensor, position your mouse pointer along the y
% (vertical) axis and click to enter 50 data points ('maxx') as you move
% the mouse pointer up and down the y axis simulating a signal. If the
% number of data points exceeds 20 ('maxdisplay'), the x-axis maximum is
% rescaled to twice that (line 32). If the data amplitude equals or exceeds
% ('maxy'), the y axis is rescaled to 1.1 times the data amplitude (line
% 36). (If you wish to change 'maxx', 'maxdisplay' or 'maxy', you may
% do so in lines 18, 19, and 20, respectively).
% 
% To do this with real time data from a sensor, replace lines 25 and 26
% with the code that acquires one data point from your sensor.
% Tom O'Haver (toh@umd.edu) 2018
clf
maxx=50; %  maximum number of data points acquired before stopping.
maxdisplay=20; % maximum number of data points before re-scaling x axis.
maxy=10; % maximum value of y axis
axis([0 maxdisplay 0 maxy]);
hold on
y=zeros(1,maxx); % y is vector of simulated data points
for n=1:maxx-1
    [clickX,clickY] = ginput(1); % Wait for a mouse click
    y(n)=clickY; % simulate a data point by the y position of the click on the graph
    if n==1
        plot(1,y(1)) % Plot the first simulated data point
        xlabel('x');ylabel('y') % Label axes
        title('Real Time Rescaling Plot Demo')
    else
        if n>maxdisplay
            maxdisplay=2*maxdisplay;
        end
        axis([0 maxdisplay 0 maxy]);  
        if y(n)>=maxy
            maxy=1.1.*y(n);
            axis([0 maxdisplay 0 maxy]);
        end
        plot([n-1 n],[y(n-1) y(n)],'k') % Draw black line from previous point
    end
end
hold off

