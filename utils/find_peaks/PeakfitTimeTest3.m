% Effect of number of points on peakfit execution time
warning off all
PeakShape=1;
NumPeaks=1;
Noise=.01;
trial=1;
clear x y time points
for interval=.005:.01:.2,
  x=[0:interval:10]';
  y=exp(-(x-5).^2)+Noise*randn(size(x));
  tic
  [FitResults,Error]=peakfit([x y],0,0,NumPeaks,PeakShape,33,1,0,0,60,0);
  et=toc;
  NumPoints=length(x);
  time(trial)=et;
  points(trial)=NumPoints;
  trial=trial+1;
  % disp(['Elapsed time ' num2str(toc)  '   NumPeaks: ' num2str(NumPeaks)   '   Error: ' num2str(Error)])
end
clf
plotit(points,time,1)
xlabel('Number of points')
ylabel('Execution time, seconds')
axis([1 2000 0 max(time)])