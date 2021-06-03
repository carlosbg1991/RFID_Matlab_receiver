function [M,A]=autopeaksplot(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Peak detection and height and area measurement for peaks of arbitrary
% shape in x,y time series data. Returns the peak table in the matrix P and
% the peak detection parameters in the vector A. Plots the data and numbers
% the peaks in Figure window 1 and plots the individual peaks in Figure
% window 2.
% M=autopeaksplot(x,y) works well in some cases, but if not try 
% M=autopeaksplot(x,y,n), using different values of n (roughly the number of
% peaks that would fit into the signal record) until it detects the peaks
% that you want to measure. For the most precise control over peak
% detection, you can specify all the peak detection parameters by typing
% M=autopeaksplot(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup)
%
% The syntax is similar syntax to autofindpeaks.m, in that that the peak
% detection parameters SlopeThreshold, AmpThreshold, smoothwidth,
% peakgroup, smoothtype can be omitted and the function will calculate
% initial values based on the number of data points and/or the optional
% peak density n (the last of 3 arguments after x and y). If you do specify
% all the peak detection parameters, they can be scalars or vectors with
% one entry for each segment (like findpeaskSG.m).
%
% Locates and measures the positive peaks in a noisy x-y time series data.
% Detects peaks by looking for downward zero-crossings in the first
% derivative whose upward slopes exceed SlopeThreshold. Returns matrix (M)
% containing peak number, the absolute peak height, peak-valley difference,
% perpendicular drop area, and the tangent skim area of each peak.
% Arguments "SlopeThreshold", "AmpThreshold" and "smoothwidth" control peak
% sensitivity of each segment. Higher values will neglect smaller features.
% "smoothwidth" is a vector of the widths of the smooths applied before
% peak detection; larger values ignore narrow peaks. If smoothwidth=0, no
% smoothing is performed. "peakgroup" is a vector of the number of points
% around the top part of the peak that are taken for initial estimate of
% the peak center and width (minimum 3). Note: this function uses
% smoothing only for peak detection; it performs measurements on the raw
% unsmoothed y data. If the data are noisy, it may be beneficial to smooth
% the y data yourself before calling autopeaks.m, using any smooth function
% of your choice.
%
% See https://terpconnect.umd.edu/~toh/spectrum/Integration.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver, 2017.  Version 1, January, 2017
%
% The script testautopeaks.m will run all of the following examples with
% a 1-second pause between each:
%
% Example 1: sin(x).^2 has theoretical peaks at x=0.5pi, 1.5pi, 2.5pi,
% 3.5pi..., with peak heights of 1.0 and peak areas of pi/2 = 1.5708.
%  disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
%  x=[0:.01:20]';y=sin(x).^2;autopeaksplot(x,y,0,0,5,5,1)
%
% Example 2: The built-in "humps" function has two peaks, no noise.
%  autopeaksplot(0:.01:2,humps(0:.01:2),0,0,1,1,1)
%
% Example 3: Series of peaks that get progressively taller and wider.
%  x=[0:.01:5]';y=x.*sin(x.^2).^2;autopeaksplot(x,y,0,0,5,5,1)
%
% Example 4: Like example 3, with random white noise added.
%  x=[0:.01:5]';y=.1.*randn(size(x))+x.*sin(-x.^2).^2;
%  autopeaksplot(x,y,.001,.5,15,15,1)
%
% Example 5; Like example 3, with added rising baseline.
%  x=[0:.01:5]';y=x+x.*sin(x.^2).^2;autopeaksplot(x,y,0,0,5,5,1)
%
% Example 6: Gaussian on linear baseline, theoretical area = 1.7725
%  x=[0:.1:10];y=2+x/10+exp(-(x-5).^2)+.01.*randn(size(x));
%  disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
%  autopeaksplot(x,y,0.0001,0.3,3,5,1)
%  Peak-valley and Tan skim are most accurate height and area measures.
%
% Example 7: Two overlapping Gaussians, zero baseline, theoretical area = 1.7725
%  x=[0:.05:10];y=exp(-(x-6).^2)+exp(-(x-3.5).^2)+.01.*randn(size(x));
%  disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
%  autopeaksplot(x,y,0.0001,0.5,8,6,1)
%  PeakMax and Perp drop are most accurate height and area measures.
%
% Example 8: Narrow Gaussian peak on sloping linear baseline, followed by
% much broader peak, theoretical heights 1 and 1, areas 1.59 and 31.9.
%  x=1:.2:150;
%  y=(150-x)./200+gaussian(x,30,1.5)+gaussian(x,90,30);
%  autopeaksplot(x,y,[0.001 .00005],[.6 .6],[5 120],[8 150],1)

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
switch nargin
    % 'nargin' is the number of input arguments
    case 1  % One argument only
        % Assumne that the argument must be a matrix of data.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(x);
        if datasize(1)<datasize(2),x=x';end
        datasize=size(x);
        if datasize(2)==1, %  Must be ipeak(y-vector)
            y=x;
            x=[1:length(x)]'; % Create an independent variable vector
        else
            % Must be ipeak(DataMatrix)
            y=x(:,2);
            x=x(:,1); % Split matrix argument
        end
        % Calculate default values of peak detection parameters
        PeakDensity=20;
        % Estimate approx number of points in a peak half-width
        WidthPoints=length(y)/PeakDensity;
        SlopeThreshold=WidthPoints^-2;
        AmpThreshold=abs(min(y)+0.5*(max(y)-min(y)));
        smoothwidth=10;
        peakgroup=10;
        if peakgroup>100,peakgroup=100;end   % Keep peakgroup below 100
        if smoothwidth>100,smoothwidth=100;end   % Keep smoothwidth below 100
        smoothtype=3;
        A=[SlopeThreshold AmpThreshold smoothwidth peakgroup];
    case 2
        PeakDensity=20;
        % Estimate approx number of points in a peak half-width
        WidthPoints=length(y)/PeakDensity;
        SlopeThreshold=WidthPoints^-2;
        AmpThreshold=abs(min(y)+0.1*(max(y)-min(y)));
        smoothwidth=round(WidthPoints/3);
        peakgroup=round(WidthPoints/3);
        if peakgroup>100,peakgroup=100;end   % Keep peakgroup below 100
        if smoothwidth>100,smoothwidth=100;end   % Keep smoothwidth below 100
        smoothtype=3;
        A=[SlopeThreshold AmpThreshold smoothwidth peakgroup];
    case 3
        % Must be separate x and y data vectors plus a peak density
        % estimate.
        % Calculate values of peak detection parameters
        % arguments based on the peak density, PeakD
        PeakDensity=SlopeThreshold;
        % Estimate approx number of points in a peak half-width
        WidthPoints=length(y)/PeakDensity;
        SlopeThreshold=WidthPoints^-2;
        AmpThreshold=abs(min(y)+0.02*(max(y)-min(y)));
        smoothwidth=round(WidthPoints/3);
        peakgroup=round(WidthPoints/3);
        if peakgroup>100,peakgroup=100;end   % Keep peakgroup below 100
        if smoothwidth>100,smoothwidth=100;end   % Keep smoothwidth below 100
        smoothtype=3;
        A=[SlopeThreshold AmpThreshold smoothwidth peakgroup];
end

if peakgroup<3,peakgroup=3;end
P=findpeaksSG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);
sizeP=size(P);
NumPeaks=sizeP(1);
clear PeakHeight
clear M
hold off
% M=zeros(NumPeaks,6);
goodpeak=1;
figure(2)
clf
for PeakNumber=1:NumPeaks,
    PeakPos=P(PeakNumber,2);
    PeakWidth=P(PeakNumber,4);
    PeakIndex=val2ind(x,PeakPos);
    x1=val2ind(x,PeakPos-1.3.*PeakWidth); % Left-most end of peak
    x2=val2ind(x,PeakPos+1.3.*PeakWidth);  % Right-most end of peak
    Halfway=val2ind(x,PeakPos); % x-value of center of peak
    miny=min(y(x1:x2)); % Minimum y value of that region
    yy(x1:x2)=y(x1:x2)-miny;  % Subtract minimum from that region
    % Calculate height, peak-valley, and perpendicular drop area
    PeakMax=y(PeakIndex);
    LeftSideMinimum=min(yy(x1:Halfway));
    RightSideMinimum=min(yy(Halfway:x2));
    LeftValleyIndex=x1+val2ind(yy(x1:Halfway),LeftSideMinimum);
    RightValleyIndex=PeakIndex+val2ind(yy(Halfway:x2),RightSideMinimum);
    if RightValleyIndex>length(x),RightValleyIndex=length(x);end
    PDarea=sum(y(LeftValleyIndex:RightValleyIndex)).*(x(2)-x(1)); % Perpendicular Drop area
    PVHeight=PeakMax-(interp1([x(LeftValleyIndex) x(RightValleyIndex)],[y(LeftValleyIndex) y(RightValleyIndex)],PeakPos));
    interpy=interp1([x(LeftValleyIndex) x(RightValleyIndex)],[y(LeftValleyIndex) y(RightValleyIndex)],x(LeftValleyIndex:RightValleyIndex));
    TSarea=sum((y(LeftValleyIndex:RightValleyIndex)-interpy).*(x(x1+1)-x(x1)));
    % [PeakNumber PeakPos PeakMax PVHeight]
    if PVHeight>AmpThreshold,
        M(goodpeak,:)=[goodpeak PeakPos  PeakMax PVHeight PDarea TSarea];
        % Plot individual peaks in a square subplot grid
        side=round(sqrt(NumPeaks)+.4); % Number of subplots on a side of square array
        if PeakNumber>side^2,
        else
            subplot(side,side,PeakNumber) % Set up side^2 square set of subplots
        end
        plot(x(x1:x2),y(x1:x2)) % Plot the individual peak
        title(['Peak ' num2str(goodpeak) ]) % Title the individual peak
        xlabel('X');ylabel('Y')
        % Mark valleys with magenta lines
        hold on
        if LeftValleyIndex<1,LeftValleyIndex=1;end
        plot([x(LeftValleyIndex) x(LeftValleyIndex)],[miny PeakMax],'m.-')
        if RightValleyIndex>length(x),RightValleyIndex=length(x);end
        plot([x(RightValleyIndex) x(RightValleyIndex)],[miny PeakMax],'m.-')
        plot([x(LeftValleyIndex) x(RightValleyIndex)],[interpy(1) interpy(length(interpy))],'c.-')
        plot(x(PeakIndex),PeakMax,'ro')
        hold off
        goodpeak=goodpeak+1;
    end
end % for PeakNumber
figure(1)
if goodpeak-1>0,
    figure(1)
    clf
    plot(x,y);
    text(M(:,2),M(:,3),num2str(M(:,1))),
    xlabel('X');
    ylabel('Y')
end

function P=findpeaksSG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Segmented peak finder, same syntax as findpeaksG except the 3rd to 6th
% input arguments can be vectors with one entry for each segment.
% function P=findpeaksSG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Locates and measures the positive peaks in a noisy x-y time series data.
% Detects peaks by looking for downward zero-crossings in the first
% derivative whose upward slopes exceed SlopeThreshold. Returns list (P)
% containing peak number and position, height, width, and area of each
% peak. Arguments "slopeThreshold", "ampThreshold" and "smoothwidth"
% control peak sensitivity of each segment. Higher values will neglect
% smaller features. "Smoothwidth" is a vector of the widths of the smooths
% applied before peak detection; larger values ignore narrow peaks. If
% smoothwidth=0, no smoothing is performed. "Peakgroup" is a vector of the
% number points around the top part of the peak that are taken for
% measurement. If Peakgroup=0 the local maximum is taken as the peak height
% and position. The argument "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar) If
%   smoothtype=2, triangular (2 passes of sliding-average) If smoothtype=3,
%   pseudo-Gaussian (3 passes of sliding-average)
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver, 2016.  Version 1, November, 2016
%
% Example: Find, measure, and plot 2 noisy peaks with very different widths
%   x=1:.2:100;
%   y=gaussian(x,20,1.5)+gaussian(x,80,30)+.02.*randn(size(x));
%   plot(x,y,'c.')
%   P=findpeaksSG(x,y,[0.001 .0001],[.2 .2],[5 10],[10 50],3)
%   text(P(:,2),P(:,3),num2str(P(:,1)))
%   disp(' ')
%   disp('           peak #    Position      Height     Width      Area')
%   disp(P)
%
% Related functions:
% findpeaksG.m, findvalleys.m, findpeaksL.m, findpeaksb.m, findpeaksb3.m,
% findpeaksplot.m, peakstats.m, findpeakSNR.m, findpeaksGSS.m,
% findpeaksLSS.m, findpeaksfit.m, findsteps.m, findsquarepulse.m, idpeaks.m

% Copyright (c) 2016 Thomas C. O'Haver
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
if nargin~=7;smoothtype=1;end  % smoothtype=1 if not specified in argument
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end
if smoothwidth<1;smoothwidth=1;end
if isscalar(AmpThreshold),AmpThreshold=AmpThreshold.*ones(size(SlopeThreshold));end
if isscalar(smoothwidth),smoothwidth=smoothwidth.*ones(size(SlopeThreshold));end
if isscalar(peakgroup),peakgroup=peakgroup.*ones(size(SlopeThreshold));end
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=SegmentedSmooth(deriv(y),smoothwidth,smoothtype);
else
    d=deriv(y);
end

P=[0 0 0 0 0];
vectorlength=length(y);
NumSegs=length(SlopeThreshold);
peak=1;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth-1,
    Seg=fix(1+NumSegs./(vectorlength./j));
    n=round(peakgroup(Seg)/2+1);
    % [j Seg]
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold(Seg), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpThreshold(Seg),  % if height of peak is larger than AmpThreshold
                xx=zeros(size(peakgroup(Seg)));yy=zeros(size(peakgroup(Seg)));
                for k=1:peakgroup(Seg), % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    % groupindex=groupindex
                    xx(k)=x(groupindex);
                    yy(k)=y(groupindex);
                end
                if peakgroup(Seg)>2,
                    [Height,Position,Width]=gaussfit(xx,yy);
                    PeakX=real(Position);   % Compute peak position and height of fitted parabola
                    PeakY=real(Height);
                    MeasuredWidth=real(Width);
                    % if the peak is too narrow for least-squares technique to work
                    % well, just use the max value of y in the sub-group of points near peak.
                else
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end
                % Construct matrix P. One row for each peak detected,
                % containing the peak number, peak position (x-value) and
                % peak height (y-value). If peak measurement fails and
                % results in NaN, or if the measured peak height is less
                % than AmpThreshold, skip this peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold(Seg),
                    % Skip this peak
                else % Otherwise count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth];
                    peak=peak+1; % Move on to next peak
                end
            end
        end
    end
end
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function d=deriv(a)
% First derivative of vector using 2-point central difference.
%  T. C. O'Haver, 1988.
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
    d(j)=(a(j+1)-a(j-1)) ./ 2;
end

function SmoothedSignal=SegmentedSmooth(y,smoothwidths,type,ends)
%   SegmentedSmooth(y,w,type,ends) divides y into a number of equal-length
%   segments defined by the length of the vector 'smoothwidths', then
%   smooths each segment with smooth of type 'type' and width defined by
%   the elements of vector 'smoothwidths'.
%    Version 1, November 2016.
%   The argument "type" determines the smooth type:
%     If type=1, rectangular (sliding-average or boxcar)
%     If type=2, triangular (2 passes of sliding-average)
%     If type=3, pseudo-Gaussian (3 passes of sliding-average)
%     If type=4, pseudo-Gaussian (4 passes of same sliding-average)
%     If type=5, multiple-width (4 passes of different sliding-averages)
%   The argument "ends" controls how the "ends" of the signal
%   (the first w/2 points and the last w/2 points) are handled.
%     If ends=0, the ends are zero. (In this mode the elapsed
%       time is the fastest and is independent of the smooth width).
%     If ends=1, the ends are smoothed with progressively
%       smaller smooths the closer to the end. (In this mode the
%       elapsed time increases with increasing smooth widths).
%   SegmentedSmooth(Y,w,type) smooths with ends=0.
%   SegmentedSmooth(Y,w) smooths with type=1 and ends=0.
%
%   Examples: 3-segment smooth of random white noise, smooth widths of
%   2,20, and 200.
%     x=1:10000;y=randn(size(x));
%     plot(x,SegmentedSmooth(y,[2 20 200],3,0))
%
%     20-segment smooth, odd smooth widths from 1 to 41:
%     plot(x,SegmentedSmooth(y,[1:2:41],3,0))

%   Copyright (c) 2012, Thomas C. O'Haver

if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
ly=length(y);
NumSegments=length(smoothwidths);
SegLength=round(ly./NumSegments);
SmoothSegment=zeros(ly,NumSegments);
SmoothedSignal=zeros(1,ly);
for Segment=1:NumSegments,
    SmoothSegment(:,Segment)=fastsmooth(y,smoothwidths(Segment),type,ends);
    startindex=(1+(Segment-1)*SegLength);
    endindix=startindex+SegLength-1;
    if endindix>ly,endindix=ly;end
    indexrange=startindex:endindix;
    SmoothedSignal(indexrange)=SmoothSegment(indexrange,Segment);
end

function SmoothY=fastsmooth(Y,w,type,ends)
% fastsmooth(Y,w,type,ends) smooths vector Y with smooth
%  of width w. Version 3.0, October 2016.
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar)
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
%   If type=4, pseudo-Gaussian (4 passes of same sliding-average)
%   If type=5, multiple-width (4 passes of different sliding-average)
% The argument "ends" controls how the "ends" of the signal
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively
%     smaller smooths the closer to the end. (In this mode the
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Examples:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
%
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%
% x=1:100;
% y=randn(size(x));
% plot(x,y,x,fastsmooth(y,5,3,1),'r')
% xlabel('Blue: white noise.    Red: smoothed white noise.')
%
% Copyright (c) 2012, Thomas C. O'Haver
%

switch type
    case 1
        SmoothY=sa(Y,w,ends);
    case 2
        SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
        SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
    case 4
        SmoothY=sa(sa(sa(sa(Y,w,ends),w,ends),w,ends),w,ends);
    case 5
        SmoothY=sa(sa(sa(sa(Y,round(1.6*w),ends),round(1.4*w),ends),round(1.2*w),ends),w,ends);
end

function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
    s(k+halfw-1)=SumPoints;
    SumPoints=SumPoints-Y(k);
    SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
        SmoothY(k)=mean(Y(1:(2*k-1)));
        SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
end

function [Height, Position, Width]=gaussfit(x,y)
% Converts y-axis to a log scale, fits a parabola
% (quadratic) to the (x,ln(y)) data, then calculates
% the position, width, and height of the
% Gaussian from the three coefficients of the
% quadratic fit.  This is accurate only if the data have
% no baseline offset (that is, trends to zero far off the
% peak) and if there are no zeros or negative values in y.
%
% Example 1: Simplest Gaussian data set
% [Height, Position, Width]=gaussfit([1 2 3],[1 2 1])
%    returns Height = 2, Position = 2, Width = 2
%
% Example 2: best fit to synthetic noisy Gaussian
% x=50:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
% [Height,Position,Width]=gaussfit(x,y)
%   returns [Height,Position,Width] clustered around 100,100,100.
%
% Example 3: plots data set as points and best-fit Gaussian as line
% x=[1 2 3 4 5];y=[1 2 2.5 2 1];
% [Height,Position,Width]=gaussfit(x,y);
% plot(x,y,'o',linspace(0,8),Height.*gaussian(linspace(0,8),Position,Width))
% Copyright (c) 2012, Thomas C. O'Haver
% To prevent problems from taking the log of zero or negative values,
% make the lowest value of y equal to 1% of the maximum value.
maxy=max(y);
for p=1:length(y),
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
logyyy=log(abs(y));
[coef,S,MU]=polyfit(x,logyyy,2);
c1=coef(3);c2=coef(2);c3=coef(1);
% Compute peak position and height or fitted parabola
Position=-((MU(2).*c2/(2*c3))-MU(1));
Height=exp(c1-c3*(c2/(2*c3))^2);
Width=norm(MU(2).*2.35703/(sqrt(2)*sqrt(-1*c3)));
