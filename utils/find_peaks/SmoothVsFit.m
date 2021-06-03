% Script that compares 3 methods for measuring the peak height of a very
% noisy peak whose peak height = 1.00 and whos signal-to-noise ratio is 1.  
x=0:.01:10;
NumTrials=100;
SmoothWidth=50;
maxy=zeros(size(NumTrials));
SmoothMaxy=zeros(size(NumTrials));
PeakfitHeight=zeros(size(NumTrials));
for Trial=1: NumTrials,
    y=exp(-(x-5).^2)+randn(size(x));
    maxy(Trial)=max(y);
    SmoothMaxy(Trial)=max(fastsmooth(y,SmoothWidth,3));
    r=peakfit([x;y]);
    PeakfitHeight(Trial)=r(3);
end
AverageMaxy=mean(maxy);
AverageSmoothMaxy=mean(SmoothMaxy);
AveragePeakfitHeight=mean(PeakfitHeight);
SDmaxy=std(maxy);
SDSmoothMaxy=std(SmoothMaxy);
SDPeakfitHeight=std(PeakfitHeight);
disp('  ')
disp(['NumTrials = ' num2str(NumTrials)  '      SmoothWidth = ' num2str(SmoothWidth)]);
disp(['      Method         Maximum y    Max Smoothed y    Peakfit '])
disp(['Average peak height   '  num2str(AverageMaxy) '        '  num2str(AverageSmoothMaxy) '         ' num2str(AveragePeakfitHeight)])
disp(['Standard deviation    '  num2str(SDmaxy) '       '  num2str(SDSmoothMaxy) '         ' num2str(SDPeakfitHeight)])