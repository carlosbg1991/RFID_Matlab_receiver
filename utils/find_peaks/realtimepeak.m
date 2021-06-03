% Here's an example of real-time peak detection that you can run on your
% computer in Matlab or Octave without any additional hardware.  It uses
% the mouse pointer to simulate an analog-to-digital converter. When you
% run this script, it will display a graph. To simulate a signal coming
% from an instrument or sensor, position your mouse pointer along the y
% (vertical) axis and click to enter 20 data points as you move the mouse
% pointer up and down the y axis simulating one or more peaks. Each time
% your mouse clicks form a peak (that is, go up and then down again), the
% program will register and label the peak on the graph and print out its x
% and y values. (If you wish you can change 'maxx', the maximum number of
% data points, in line 17).
% 
% To do this with real time data from a sensor, replace lines 22 and 23
% with the code that acquires one data point from your sensor.
% Tom O'Haver (toh@umd.edu) 2018
clf
maxx=20; % Sets the maxumum number of data points
axis([0 20 0 10]);
hold on
y=zeros(1,maxx); % y is vector of simulated data points
for n=1:maxx-1
    [clickX,clickY] = ginput(1); % Wait for a mouse click
    y(n)=clickY; % simulate a data point by the y position of the click on the graph
    if n==1
        plot(1,y(1)) % Plot the first simulated data point
        xlabel('x');ylabel('y') % Label axes
        title('Real Time Peak Detection Demo')
    else
        plot([n-1 n],[y(n-1) y(n)],'k') % Draw black line from previous point
        if n>2 % Start peak detection when 3 points have been acquired.
            if y(n-1)>y(n-2) % If a point is greater that the the previous one
                if y(n-1)>y(n) % AND greater than the following one, register a peak.
                    disp(['Peak detected at x=' num2str(n-1) ' and y=' num2str(y(n-1))])
                    text(n-1,y(n-1),' peak') % Label the peak on the graph
                end
            end
        end
    end
end
hold off

