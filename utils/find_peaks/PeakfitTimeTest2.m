% Effect of peak shape on peakfit execution time
load DataMatrix3
warning off all
NumPeaks=2; % Allowed values 1-4
NumTrials=10;
disp(' ')
disp('Comparison of Gaussian peak shape variations in peakfit function:')
disp(['Number of peaks = ' num2str(NumPeaks)])
disp(' ')
PeakShape=1;
tic;[FitResults,Error]=peakfit(DataMatrix3,1855.5,387,NumPeaks,PeakShape,33,NumTrials,0,0,60,1);et(1)=toc;er(1)=Error(1);
disp('Gaussian')
disp(['Elapsed time ' num2str(et(1))  '   Error: ' num2str(Error)])
drawnow
disp(' ')
PeakShape=6;
tic;[FitResults,Error]=peakfit(DataMatrix3,1855.5,387,NumPeaks,PeakShape,33,NumTrials,0,0,60,1);et(2)=toc;er(2)=Error(1);
disp('Equal-width Gaussians.')
disp(['Elapsed time ' num2str(et(2))  '   Error: ' num2str(Error)])
drawnow
disp(' ')
PeakShape=11;
tic;[FitResults,Error]=peakfit(DataMatrix3,1855.5,387,NumPeaks,PeakShape,33,NumTrials,0,0,60.*ones(1,NumPeaks),1);et(3)=toc;er(3)=Error(1);
disp('Fixed-width Gaussians.')
disp(['Elapsed time ' num2str(et(3))  '   Error: ' num2str(Error)])
drawnow
disp(' ')
PeakShape=5;
tic;[FitResults,Error]=peakfit(DataMatrix3,1855.5,387,NumPeaks,PeakShape,33,NumTrials,0,0,60,1);et(4)=toc;er(4)=Error(1);
disp('Exponentially broadened Gaussians with default start.')
disp(['Elapsed time ' num2str(et(4))  '   Error: ' num2str(Error)])
drawnow
disp(' ')
PeakShape=8;
tic;[FitResults,Error]=peakfit(DataMatrix3,1855.5,387,NumPeaks,PeakShape,33,NumTrials,0,0,60,1);et(5)=toc;er(5)=Error(1);
disp('Exponentially broadened equal-width Gaussians with default start.')
disp(['Elapsed time ' num2str(et(5))  '   Error: ' num2str(er(5))])
drawnow
disp(' ')
PeakShape=5;
if NumPeaks==1,start=[1814 66];end
if NumPeaks==2,start= [1814 66 1920 77];end
if NumPeaks==3,start=[1814 66 1920 77 2400.5 59];end
if NumPeaks==4,start=[1814 66 1920 77 2400.5 59 2400.5 59];end
tic;[FitResults,MeanFitError]=peakfit(DataMatrix3,1855.5,387,NumPeaks,PeakShape,33.0226,NumTrials, start, 0, 0 );et(6)=toc;er(6)=Error(1);
disp('Exponentially broadened Gaussians with good start.')
disp(['Elapsed time ' num2str(et(6))  '   Error: ' num2str(er(6))])

