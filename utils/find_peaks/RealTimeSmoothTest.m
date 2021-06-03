% Example of real-time smoothing. Set SmoothWidth in line 18. PlottingOn=1
% to see plots. Change 'maxx' in line 15 to control the total number of
% data points and SmoothWidth in line 19). Uncomment either line 24
% (rectangular), 25 (triangular), or 26 (Gaussian), or write a smooth
% vector of your choice. To do this with real time data from a sensor,
% comment out line 31 and replace line 35 with the code that acquires one
% data point from your sensor. For the fastest data acquisition speeds, set
% PlottingOn (line 20) to 0. Requires triangle and gaussian functions. Tom
% O'Haver (toh@umd.edu) 2019. Version 2
%
%
format compact
format short g
clf;clear
maxx=500; % Sets the maxumum number of data points
hold on
grid
noise=.02;
SmoothWidth=21; % Smooth width (odd integer)
PlottingOn=1; % 0=no plot (much faster)  1=draw plot (slower)
% SmoothWidth=makeodd(SmoothWidth); % Optional to insure oddness.
% Uncomment one of the following lines to define desired smooth shape.
SmoothVector=ones(1,SmoothWidth);SmoothVector=SmoothVector./(sum(SmoothVector));
% SmoothVector=triangle(1:SmoothWidth,SmoothWidth/2,SmoothWidth/2);SmoothVector=SmoothVector./(sum(SmoothVector));
% SmoothVector=gaussian(1:SmoothWidth,SmoothWidth/2,SmoothWidth/2);SmoothVector=SmoothVector./(sum(SmoothVector));
% SmoothVector=1-exp(-(1:SmoothWidth)./SmoothWidth.*10);SmoothVector=SmoothVector./(sum(SmoothVector));
y=zeros(1,maxx); % y is vector of simulated data points
sy=zeros(1,maxx); % sy is vector of smoothed data points
RawData=rectanglepulse(1:maxx,maxx/2,maxx/3)+noise.*randn(size(y)); % Generate rectangular sumulated data
% RawData=gaussian(1:maxx,maxx/2,maxx/3)+noise.*randn(size(y)); % Generate Gaussian sumulated data
axis([0 maxx -(.1+3*noise) 1.1+3*noise]); % Change to suit your data

tic
for n=1:maxx % Stop after maxx data ppints
    y(n)=RawData(n);  % <<< replace this line with call to your data acquisition harware
    FirstPoint=n-SmoothWidth;
    if FirstPoint<1;FirstPoint=1;end
    LastPoint=n-1;
    if LastPoint<1;LastPoint=1;end
    if n>SmoothWidth
        SmoothedPoint=sum(SmoothVector.*(y(FirstPoint:LastPoint)));
        sy(n)=SmoothedPoint;
        % Add statements here to output smoothed data "SmoothedPoint".
    end
    if PlottingOn==1
        if n==1 %
            plt = Plot(1,sy(1)); % Plot the first simulated data point
            xlabel(['Time      SmoothWidth= ' num2str(SmoothWidth)]);
            ylabel('Signal amplitude')
            title('Real Time Smooth Demo.  Blue=original   Red=Smoothed')
        else
            LineStart=n-1;
            if LineStart<1;LineStart=1;end
            LineEnd=n;
            if n>SmoothWidth
                plt = Plot([n-1 n],[sy(LineStart) sy(LineEnd)],'r'); % Draw smoothed signal as red
            end
            plt = Plot([n-1 n],[y(LineStart) y(LineEnd)],'b'); % Draw raw signal as blue
        end % if n==1
        drawnow % Comment this line out to delay plotting until end of signal
    end % plot
end % for n=1:maxx
hold off
ElapsedTime=toc;
TimePerPoint=ElapsedTime./maxx;
CyclesPerSecond=1/TimePerPoint;
disp('        Plot?    NumPoints    PointsPerSecond')
disp([PlottingOn  maxx CyclesPerSecond])
if PlottingOn==0
    plot(1:n,y,1:n,sy) % Post-run plotting
    grid
    xlabel(['Time      SmoothWidth= ' num2str(SmoothWidth)]);
    ylabel('Signal amplitude')
    title('Real Time Smooth, post-run plot.  Blue=original   Red=Smoothed')
end



