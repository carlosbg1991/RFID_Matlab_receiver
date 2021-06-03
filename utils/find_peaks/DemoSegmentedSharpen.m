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
y=gaussian(x,6,2)+gaussian(x,14,3)+gaussian(x,22,4)+gaussian(x,30,5);
NumSegments=30; % Sets the number of segmemnts
start1=25;
end1=100;
factor1=(start1:(end1-start1)/NumSegments:end1).^2; % Increasing values of factor 1
start2=30;
end2=110;
factor2=(start2:(end2-start2)/NumSegments:end2).^4; % Increasing values of factor 2
smoothwidth=1.*ones(size(factor1)); % Same smooth width for all peaks
ey=SegmentedSharpen(y,factor1,factor2,smoothwidth);
plot(x,y,x,ey,'r')
axis([ 0 40 -.1 1.5])
title('Demonstration of segmented sharpen function.')
ylabel('y')
xlabel('Four Gaussian peaks (blue) with increasing widths subjected to increasing sharpening (red)')

OriginalWidths=[halfwidth(x,y,6)  halfwidth(x,y,14)  halfwidth(x,y,22)  halfwidth(x,y,30)];
SharpenedWidths=[halfwidth(x,ey,6)  halfwidth(x,ey,14)  halfwidth(x,ey,22)  halfwidth(x,ey,30)];
PercentSharpened=100.*(OriginalWidths-SharpenedWidths)./OriginalWidths;
disp(' ')
disp('Effect on peak width')
disp('                   Peak 1      Peak 2       Peak 3       Peak 4')
disp(['Original Widths:     '  num2str(OriginalWidths) ])
disp(['Sharpened Widths:    '  num2str(SharpenedWidths) ])
disp(['Percent Sharpened:   '  num2str(PercentSharpened) ])
disp(' ')
disp('Percent change in individual peak areas measured by perpendicular drop:')
PercentChangeInTotalArea=100*(sum(y)-sum(ey))/sum(y);
Py=autopeaks(x,y,5);
Pey=autopeaks(x,ey,5);
PercentDifference=100.*(Pey-Py)./Py;
disp('       Peak 1      Peak 2       Peak 3       Peak 4')
disp(PercentDifference(:,5)')
disp(['Change In Total Area = ' num2str(PercentChangeInTotalArea) '%' ])