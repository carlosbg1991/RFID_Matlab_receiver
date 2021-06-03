% Here's an example of real-time smoothing that you can run on your
% computer in Matlab or Octave without any additional hardware.  Line 21
% generates a simulated signal and a smoothing vector. Lines 23-37 perform
% the actual point-by-point smoothing operation. (If you wish you can
% change 'maxx' in line 16 to control the total number of data points and
% SmoothWidth in line 17). Set the smooth shape in line 18 (1=sliding
% average, 2=trialngular, 3=Gaussian). To do this with real time data from a
% sensor, comment out line 20 and replace line 24 with the code that
% acquires one data point from your sensor.
% Requires triangle.m, gaussian.m, and fastsmooth.m functions. Tom O'Haver
% (toh@umd.edu) 2019
%
format compact
format short g
clf;clear
maxx=300; % Sets the maxumum number of data points
SmoothWidth=41; % Smooth width (odd integer)
SmoothType=1; % 1=sliding average, 2=triangular, 3=Gaussian
noise=.5; % standard deviation of normally distributed white noise
PlottingOn=1; % 0=no plot (much faster)  1=draw plot (slower)
RawData=rectanglepulse(1:maxx,maxx/2,maxx/3)+noise.*randn(1,maxx); % Generate simulated data
tic
TempSmoothWidth=1;
for n=1:maxx
    y(n)=RawData(n);  % <<< replace this line with call to your data acquisition harware
    if n==1
        sy(1)=0;
    else
        sy=fastsmooth(y,TempSmoothWidth,SmoothType,1);
        TempSmoothWidth=TempSmoothWidth+1;
        if TempSmoothWidth>SmoothWidth,TempSmoothWidth=SmoothWidth;end
        if PlottingOn==1
            plot(1:n,y,'.b',1:n,sy,'.r')
            axis([0 maxx min(RawData) max(RawData)])
            grid
            drawnow % Comment this line out to delay plotting until end of signal
        end
    end
end
ElapsedTime=toc;
xlabel(['Time      SmoothWidth= ' num2str(SmoothWidth)]);
ylabel('Response') % Label axes
title('Real Time Smooth 2 Demo.  Blue=original   Red=Smoothed')
% OnTime=val2ind(sy(1:250),.5)
% OffTime=250+val2ind(sy(250:500),.5)
TimePerPoint=ElapsedTime./maxx
CyclesPerSecond=1/TimePerPoint


