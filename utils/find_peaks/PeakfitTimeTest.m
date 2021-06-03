% Effect of number of peaks on peakfit execution time
load DataMatrix2
warning off all
window=450;
PeakShape=1;
peakfit(DataMatrix2,2425,window,3,PeakShape,33,1,0,0,60,1);
drawnow
for trial=1:30,
   NumPeaks=round(trial/5);
   if NumPeaks==0;NumPeaks=1;end
  tic;
  [FitResults,Error]=peakfit(DataMatrix2,2425,window,NumPeaks,PeakShape,33,1,0,0,60,0);
  time(trial)=toc;
  Peaks(trial)=NumPeaks;
  disp(['Elapsed time ' num2str(toc)  '   NumPeaks: ' num2str(NumPeaks)   '   Error: ' num2str(Error)])
factor=toc/((NumPeaks.^2)*window);
end
subplot(1,2,1)
plotit(Peaks,time,2)
title('Quadratic fit to NumPeaks')
xlabel('Number of peaks')
ylabel('Execution time, seconds')
subplot(1,2,2)
plotit(Peaks.^2,time,1)
title('Linear fit to NumPeaks^2')
xlabel('Number of peaks sqsuared')
ylabel('Execution time, seconds')