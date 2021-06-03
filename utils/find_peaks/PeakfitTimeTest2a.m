% Effect of peak shape on peakfit execution time

warning off all
Noise=0.1;
NumPeaks=4;
x=[0:.01:18];
y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2)+Noise.*randn(size(x));
y1=ExpBroaden(y',-50);
DataMatrix=[x;y1'];
plot(x,y1)

PeakShape=1;
disp(' ')
disp('Comparison of Gaussian peak shape variations in peakfit function:')
tic;
[FitResults,Error]=peakfit(DataMatrix,0,0,NumPeaks,PeakShape,50,1,0,0,1.7,0);
et(1)=toc;
er(1)=Error;
disp('Gaussian')
disp(['Elapsed time ' num2str(et(1))  '   Error: ' num2str(Error)])
drawnow

PeakShape=6;
tic;
[FitResults,Error]=peakfit(DataMatrix,0,0,NumPeaks,PeakShape,50,1,0,0,1.7,0);
et(2)=toc;
er(2)=Error;
disp('Equal-width Gaussians.')
disp(['Elapsed time ' num2str(et(2))  '   Error: ' num2str(Error)])
drawnow

PeakShape=11;
tic;
[FitResults,Error]=peakfit(DataMatrix,0,0,NumPeaks,PeakShape,50,1,0,0,1.7,0);
et(3)=toc;
er(3)=Error;
disp('Fixed-width Gaussians.')
disp(['Elapsed time ' num2str(et(3))  '   Error: ' num2str(Error)])
drawnow

PeakShape=5;
tic;
[FitResults,Error]=peakfit(DataMatrix,0,0,NumPeaks,PeakShape,50,1,0,0,1.7,0);
et(4)=toc;
er(4)=Error;
disp('Exponentially broadened Gaussians with default start.')
disp(['Elapsed time ' num2str(et(4))  '   Error: ' num2str(Error)])
drawnow

PeakShape=8;
tic;
[FitResults,Error]=peakfit(DataMatrix,0,0,NumPeaks,PeakShape,50,1,0,0,1.7,0);
et(5)=toc;
er(5)=Error;
disp('Exponentially broadened equal-width Gaussians with default start.')
disp(['Elapsed time ' num2str(et(5))  '   Error: ' num2str(er(5))])
drawnow

PeakShape=5;
start=[4 1.7 9 1.7 12 1.7 14 1.7];
tic;
[FitResults,MeanFitError]=peakfit(DataMatrix,0,0,NumPeaks,PeakShape,50,1,start,0,0);
et(6)=toc;
er(6)=Error;
disp('Exponentially broadened Gaussians with good start.')
disp(['Elapsed time ' num2str(et(6))  '   Error: ' num2str(er(6))])

