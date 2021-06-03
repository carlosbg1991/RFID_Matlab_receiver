% Here's an example of real-time signal plotting that you can run on your
% computer in Matlab or Octave without any additional hardware.  It uses
% stored data, accessed point-by-point, to simulate real-time data. If the
% number of data points exceeds 'maxdisplay', the x-axis maximum is
% rescaled to 1.1. times that (line 34). If the data amplitude equals or
% exceeds 'maxy', the y axis is rescaled to 1.1 times the data amplitude
% (line 38). (If you wish to change 'maxn', 'maxdisplay' or 'maxy', you may
% do so in lines 15-17).
% 
% Tom O'Haver (toh@umd.edu) 2018
figure(1)
clf
load DataMatrix2
maxn=1200; %  maximum number of data points acquired before stopping.
maxdisplay=1000; % maximum number of data points before re-scaling x axis.
maxy=DataMatrix2(1,2); % maximum value of y axis
miny=-.05;
minx=min(DataMatrix2(1:maxn,1));
maxx=max(DataMatrix2(1:maxn,1));
axis([0 maxx miny maxy]);
hold on
x=zeros(1,maxn); % x is vector of simulated data points x=DataMatrix2(:,1);
y=zeros(1,maxn); % y is vector of simulated data points y=DataMatrix2(:,2);
tic
for n=1:maxn-1
    x(n)=DataMatrix2(n,1);
    y(n)=DataMatrix2(n,2); % simulate a single data point from the data source
    if n==1
        plot(x(1),y(1)) % Plot the first simulated data point
        xlabel('x');ylabel('y') % Label axes
        title('Real Time Rescaling Plot Demo')
    else
        if n>maxdisplay
            maxdisplay=1.1*maxdisplay;
            axis([0 maxx miny maxy]);  
        end
        if y(n)>=maxy
            maxy=1.1.*y(n);
            axis([0 maxx miny maxy]);
        end
        plot(x(n),y(n),'.k') % Draw black dot
    end
    drawnow
end
hold off
elapsedtime=toc
TimePerPoint=elapsedtime/maxn
