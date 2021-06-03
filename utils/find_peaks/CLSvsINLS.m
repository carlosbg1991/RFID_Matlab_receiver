% Iterative non-linear least-squares (INLS) compared to Classical Least
% Squares (CLS) for the measurement of peak heights of overlapping peaks
% (Requires cls.m and peakfit.m in the Matlab path)
%
%  Create a noisy model signal
format short g
format compact
warning off all
x=200:.1:800;
NumPeaks=3;
PeakShape=1;
PeakHeights=[5 30 20];
PeakPositions=[300 500 700];
PeakWidths=[100 100 100];
extra=0;
model=modelpeaks(x,NumPeaks,PeakShape,PeakHeights,PeakPositions,PeakWidths,extra);
y=model+randn(size(x));
disp('------------------------------------------------')
% Use the peakfit.m  function to measure the peak heights in the noisy
% data, knowing the only NumPeaks and PeakShape (and, for peakshape 16 and
% 17, the peak positions)
% start=[300 100 500 150 700 100]; % start for 
% start=[100 150 700];
start=0;
tic
[FitResults,MeanFitError]=peakfit([x;y],500,600,NumPeaks,1,0,1,start,0,PeakPositions);
INLStims=toc;
INLSAccuracy=mean(abs(100.*(PeakHeights-FitResults(:,3)')./PeakHeights));
% Use the cls.m function to measure the peak heights in the noisy data.
% knowing the NumPeaks,PeakShape,PeakPositions,PeakWidths.
figure(2)
tic
h=cls(x,y,NumPeaks,PeakShape,PeakPositions,PeakWidths,extra);
CLStime=toc;
plot(x,y,x,model,'r');
MeasuredHeights=h;
CLSAccuracy=mean(abs(100.*(PeakHeights-MeasuredHeights)./PeakHeights));
disp(' Method    Positions   Widths        Execution time     Accuracy')
disp(['  CLS       known      known          ' num2str(CLStime) '        ' num2str(CLSAccuracy) ] )
disp(['  INLS      unknown    unknown        ' num2str(INLStims) '           ' num2str(INLSAccuracy) ] )

tic
[FitResults,MeanFitError]=peakfit([x;y],500,600,NumPeaks,16,0,1,start,0,PeakPositions);
INLStims=toc;
INLSAccuracy=mean(abs(100.*(PeakHeights-FitResults(:,3)')./PeakHeights));
disp(['  INLS      known      unknown        ' num2str(INLStims) '           ' num2str(INLSAccuracy) ] )

tic
[FitResults,MeanFitError]=peakfit([x;y],500,600,NumPeaks,6,0,1,start,0,PeakPositions);
INLStims=toc;
INLSAccuracy=mean(abs(100.*(PeakHeights-FitResults(:,3)')./PeakHeights));
disp(['  INLS      unknown    known          ' num2str(INLStims) '           ' num2str(INLSAccuracy) ] )

tic
w=mean(PeakWidths);
[FitResults,MeanFitError]=peakfit([x;y],500,600,NumPeaks,11,0,1,start,0,[w w w]);
INLStims=toc;
INLSAccuracy=mean(abs(100.*(PeakHeights-FitResults(:,3)')./PeakHeights));
disp(['  INLS      unknown    known (equal)  ' num2str(INLStims) '          ' num2str(INLSAccuracy) ] )
