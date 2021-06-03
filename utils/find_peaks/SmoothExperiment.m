x=[0:.1:10]';
y=exp(-(x-5).^2);
plot(x,y)
ysmoothed=fastsmooth(y,11,3,1);
plot(x,y,x,ysmoothed,'r')
disp('Fit to unsmoothed peak:')
[FitResults,FitError]=peakfit([x y])
disp('Fit to smoothed peak:')
[FitResults,FitError]=peakfit([x ysmoothed])