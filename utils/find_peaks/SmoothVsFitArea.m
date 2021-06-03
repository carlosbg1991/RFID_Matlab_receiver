% Script that compares 3 methods for measuring the peak AREA of a very
% noisy peak whose peak area = 1.772 and whos signal-to-noise ratio is 1. 
% Example:
% NumTrials = 100      SmoothWidth = 50
% Method               Area y    Area Smoothed y    Peakfit 
% Average peak area     1.7889        1.7806         1.7806
% Standard deviation    0.29035       0.27607         0.18283
increment=.01;
x=0:increment:10;
NumTrials=100;
SmoothWidth=50;
Areay=zeros(size(NumTrials));
SmoothAreay=zeros(size(NumTrials));
PeakfitArea=zeros(size(NumTrials));
for Trial=1: NumTrials,
    y=exp(-(x-5).^2)+randn(size(x));
    Areay(Trial)=sum(y)*increment;
    SmoothAreay(Trial)=sum(fastsmooth(y,SmoothWidth,3)*increment);
    r=peakfit([x;y]);
    PeakfitArea(Trial)=r(5);
end
AverageAreay=mean(Areay);
AverageSmoothMaxy=mean(SmoothAreay);
AveragePeakfitArea=mean(PeakfitArea);
SDAreay=std(Areay);
SDSmoothAreay=std(SmoothAreay);
SDPeakfitArea=std(PeakfitArea);
disp('  ')
disp(['NumTrials = ' num2str(NumTrials)  '      SmoothWidth = ' num2str(SmoothWidth)]);
disp(['Method               Area y    Area Smoothed y    Peakfit '])
disp(['Average peak area     '  num2str(AverageAreay) '        '  num2str(AveragePeakfitArea) '         ' num2str(AveragePeakfitArea)])
disp(['Standard deviation    '  num2str(SDAreay) '       '  num2str(SDSmoothAreay) '        ' num2str(SDPeakfitArea)])