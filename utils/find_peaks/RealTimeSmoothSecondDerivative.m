% Here's an example of real-time second differentiation that you can run on
% your computer in Matlab or Octave without any additional hardware.  Lines
% 9-25 generate a simulated signal and a smoothing vector. Lines 31-52
% perform the actual point-by-point smoothing operation. (If you wish you
% can change 'maxx' in line 17 to control the total number of data points
% and SmoothWidth in line 20). Uncomment either line 23 (rectangular), 24
%  (triangular), or 24 (Gaussian), or write a smooth vector of your choice.
% To do this with real time data from a sensor, comment out line 28 and
% replace line 31 with the code that acquires one data point from your
% sensor. For the fastest data acquisition speeds, disable point-by-point
% plotting by commenting out or cutting lines 18-19, 40-50, and 52. 
% Requires triangle and gaussian functions. Tom O'Haver (toh@umd.edu) 2018
format compact
format short g
clf
clear
maxx=500; % Sets the maxumum number of data points
axis([0 maxx -1 1]); % Change to suit yoour data
hold on
SmoothWidth=101; % Smooth width (odd integer)
% SmoothWidth=makeodd(SmoothWidth); % Optional to insure oddness.
% >>>Uncomment one of the following 3 lines to define desired smooth shape.<<
% SmoothVector=ones(1,SmoothWidth);SmoothVector=SmoothVector./(sum(SmoothVector));
% SmoothVector=triangle(1:SmoothWidth,SmoothWidth/2,SmoothWidth/2);SmoothVector=SmoothVector./(sum(SmoothVector));
SmoothVector=gaussian(1:SmoothWidth,SmoothWidth/2,SmoothWidth/2);SmoothVector=SmoothVector./(sum(SmoothVector));
y=zeros(1,maxx); % y is vector of simulated data points
sy=zeros(1,maxx); % sy is vector of smoothed data points
RawData=gaussian(1:maxx,maxx/2,maxx/5)+.01.*randn(size(y)); % Generate sumulated data
grid
tic
for n=SmoothWidth:maxx
    y(n)=RawData(n);  % <<< replace this line with call to your data acquisition harware
    FirstPoint=n-SmoothWidth;
    if FirstPoint<1;FirstPoint=1;end
    LastPoint=n-1;
    if LastPoint<1;LastPoint=1;end
    if n>SmoothWidth+1
        SmoothedPoint=sum(SmoothVector.*(y(FirstPoint:LastPoint)));
        sy(n)=SmoothedPoint;
    end
    xlabel(['X      SmoothWidth= ' num2str(SmoothWidth)]);ylabel('Y') % Label axes
    title('Real Time Second Derivative Demo.  Black=original   Red=Smoothed')
    if n<SmoothWidth % Comment out or remove this 'if' loop to disable plotting
        plot(1,sy(1)) % Plot the first simulated data point
    else
        plot([n-1 n],[y(n-1) y(n)],'k')
        dy(n)=maxx.*(-2.*sy(n-2)+sy(n-4)+sy(n));
        dy(1:SmoothWidth+3)=0; % Zeros partially smoothed derivative points
        plot([n-1 n],[dy(n-1) dy(n)],'r')
        drawnow % Comment out this line to disable point-by-point-plotting
    end
end
hold off
ElapsedTime=toc;
TimePerPoint=ElapsedTime./maxx

% On a standard desktop PC (Intel Core i5 3 Ghz) running Windows 10 home,
% the elapsed time per data point is 2 microseconds (without plotting)
% and 1.4 milliseconds per point with point-by-point plotting (lines 40-50).


