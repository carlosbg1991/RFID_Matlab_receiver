% Here's an example of real-time peak fitting that you can run on your
% computer in Matlab or Octave without any additional hardware.  It uses
% stored data, accessed point-by-point, to simulate real-time data. If the
% number of data points exceeds 'maxdisplay', the x-axis maximum is
% rescaled to 1.1. times that (line 34). If the data amplitude equals or
% exceeds 'maxy', the y axis is rescaled to 1.1 times the data amplitude
% (line 38). (If you wish to change 'maxn', 'maxdisplay' or 'maxy', you may
% do so in lines 14-17).
% 
% Tom O'Haver (toh@umd.edu) 2018
format compact
format short g
figure(1)
clf
load DataMatrix4
maxn=450; %  maximum number of data points acquired before stopping.
startn=150; % Number of data points skipped in the begining
nincrement=1;
maxdisplay=1000; % maximum number of data points before re-scaling x axis.
WindowWidth=75; % The number of data points over which each curve-fit is performed
numpeaks=1;
peakshape=1;
autozero=1;
maxy=1; % maximum value of y axis
miny=-.05;
minx=200;
maxx=max(DataMatrix4(1:maxn,1));
axis([minx maxx miny maxy]);
hold on
x=zeros(1,maxn); % x is vector of simulated data points x=DataMatrix2(:,1);
y=zeros(1,maxn); % y is vector of simulated data points y=DataMatrix2(:,2);
tic
for n=startn:maxn-1:nincrement
    x(n)=DataMatrix4(n,1);
    y(n)=DataMatrix4(n,2); % simulate a single data point from the data source
    if n==startn
        plot(x(1),y(1)) % Plot the first simulated data point
        xlabel('x');ylabel('y') % Label axes
        title('Real Time Peak Fitting Demo')
    else   
        n1=n-WindowWidth;if n1<1,n1=1;end
        % call peakfit funciton with plotting supressed
       figure(2);[FitResults,GOF]=peakfit([x(n1:n);y(n1:n)],0,0,numpeaks,peakshape,0,1,0,autozero,0,1);
       drawnow
      R2=GOF(2);
       figure(1);
       plot(x(n),y(n),'.k',x(n-round(WindowWidth/2)),R2,'.r') % Draw black dot
       
    end
    drawnow
end
hold off
elapsedtime=toc
TimePerPoint=elapsedtime/maxn
