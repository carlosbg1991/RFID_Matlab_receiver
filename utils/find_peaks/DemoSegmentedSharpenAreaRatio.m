% Demonstration of segmented sharpen function.
% Four Gaussian peaks with increasing widths subjected to sharpening using
% the SegmentedSharpen.m function. The script calculates the factor1,
% factor2, and smoothwidth vectors required by the SegmentedSharpen
% function (line 23) from the specified number of segments (line 15) and
% the start and end values of each vector (lines 16-20). See also
% DemoSegmentedSharpen2.m. 
% You must have SegmentedSharpen.m and halfwidth.m in the Matlab/Octave
% path. Version 1.2 June 2017. Tom O'Haver toh@umd.edu
format compact
format short g
figure(1)
clf
x=2:.01:38;
tau=-40;
width=3;
t1=10;
t2=15;
SeparationWidthRatio=(t2-t1)/width;
y=(expgaussian(x,t1,width,tau)+expgaussian(x,t2,width,tau))';
NumSegments=2; % Sets the number of segmemnts
factor1=[60 60].^2;
factor2=[60 60].^4;
smoothwidth=1.*ones(size(factor1)); % Same smooth width for all peaks
ey=SegmentedSharpen(y,factor1,factor2,smoothwidth);
plot(x,y,x,ey,'r')
axis([ 0 40 -.15 1.8])
title('Effect of peak sharpening on the measured area ratio of overlapping chiral peaks.')
ylabel('y')
xlabel('Two Gaussian peaks (blue) with equal height and width, subjected to sharpening (red)')

OriginalWidths=[halfwidth(x,y,t1)  halfwidth(x,y,t2)];
SharpenedWidths=[halfwidth(x,ey,t1)  halfwidth(x,ey,t2)];
PercentSharpened=100.*(OriginalWidths-SharpenedWidths)./OriginalWidths;
disp(' ')
disp('Effect on peak width')
disp('                   Peak 1      Peak 2')
disp(['Original Widths:     '  num2str(OriginalWidths) ])
disp(['Sharpened Widths:    '  num2str(SharpenedWidths) ])
disp(['Width Ratio:               '  num2str(OriginalWidths./SharpenedWidths) ])
disp(' ')
disp('Peak areas measured by perpendicular drop:')
PercentChangeInTotalArea=100*(sum(y)-sum(ey))/sum(y);
Py=autopeaks(x,y,2e-006,0.2,10,10,11)';
Pey=autopeaks(x,ey,2e-006,0.2,10,10,11)';
PerpDropAreas=Py(5,:);
PerpDropAreasSharpened=Pey(5,:);
RawRatio=PerpDropAreas(1)./PerpDropAreas(2);
SharpRatio=PerpDropAreasSharpened(1)./PerpDropAreasSharpened(2);
disp('           Peak 1      Peak 2    Peak1/Peak2 Ratio')
disp(['Raw data =  ' num2str([PerpDropAreas  RawRatio]) ])
disp(['Sharpened = ' num2str([PerpDropAreasSharpened   SharpRatio]) ])
text(20,1.45,['Peak separation:     '  num2str(t2-t1) ])
text(20,1.4,['factor1:     '  num2str(factor1) ])
text(20,1.35,['factor2:     '  num2str(factor2) ])

