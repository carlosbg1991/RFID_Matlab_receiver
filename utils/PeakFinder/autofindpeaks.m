function [P,DetectionParameters]=autofindpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% autofindpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype) 
% Automatic peak finder that locates and measures the positive peaks in a
% noisy x-y time series data. similar to findpeaksG or findpeaskSG except
% that the peak detection parameters SlopeThreshold, AmpThreshold,
% smoothwidth, peakgroup, smoothtype can be omitted and the function will
% calculate trial values based on the number of data points and/or the
% optional peak density (the last of 3 arguments after x and y). Returns a
% table (P) of peak number, position, absolute peak height, peak-valley
% difference, perpendicular drop area, and tangent skim area of each peak.
% Detects peaks by looking for downward zero-crossings in the first
% derivative whose upward slopes exceed SlopeThreshold. If Peakgroup=0 the
% local maximum is taken as the peak height and position. For best results,
% remove the background from the data before using this function. Optional
% input arguments "slopeThreshold", "ampThreshold" and "smoothwidth"
% control peak sensitivity of each segment. Higher values will neglect
% smaller features. "Smoothwidth" is a vector of the widths of the smooths
% applied before peak detection; larger values ignore narrow peaks. If
% smoothwidth=0, no smoothing is performed. "Peakgroup" is a vector of the
% number points around the top part of the peak that are taken for
% measurement. The argument "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar) If
%   smoothtype=2, triangular (2 passes of sliding-average) If smoothtype=3,
%   pseudo-Gaussian (3 passes of sliding-average)
% Run testautopeaks.m to test and demonstrate this function.
% In version 1.1, [P,DetectionParameters]=autofindpeaks...also returns the
% peak detection parameters as a 4-element row vector.
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver, 2016.  Version 1.1, February, 2017
%
% The script testautofindpeaks.m runs all the examples below, additionally
% plotting the data and numbering the peaks (like autofindpeaksplot.m)
%
% Example 1:  One input argument; data in single vector
% x=[0:.01:5];y=sin(10*x);autofindpeaks(y);
%
% Example 2:  One input argument; data in two columns of a matrix
% x=[0:.01:5]';y=x.*sin(x.^2).^2;M=[x y];autofindpeaks(M);
%
% Example 3: Two input arguments; data in separate x and y vectors
% x=[0:.1:100];y=(x.*sin(x)).^2;autofindpeaks(x,y);
% or  x=[0:.005:2];y=humps(x);P=autofindpeaks(x,y)
%
% Example 4:  Additional input argument (after the x,y data) to control
% peak sensitivity; higher numbers for more peaks:
%  x=[0:.1:10];y=5+5.*sin(x)+randn(size(x));autofindpeaks(x,y,3);
%  or x=[0:.1:100];y=5+5.*cos(x)+randn(size(x));autofindpeaks(x,y,10);
%  or x=[0:.1:1000];y=5+5.*cos(x)+randn(size(x));autofindpeaks(x,y,100)
%
% Example 5: Seven input arguments. Specify all peak detections parameters
% x=1:.2:100;
% y=gaussian(x,40,10)+gaussian(x,50,10)+.01.*randn(size(x));
% autofindpeaks(x,y,0.00026015,0.031007,19,21,3)
%
% Example 6: Seven input arguments. Specify all peak detections parameters, in
% this case using vectors to optimize for peaks with very different widths.
%   x=1:.2:100;
%   y=gaussian(x,20,1.5)+gaussian(x,80,30)+.02.*randn(size(x));
%   plot(x,y,'c.')
%   P=autofindpeaks(x,y,[0.001 .0001],[.2 .2],[5 10],[10 100],3)
%   text(P(:,2),P(:,3),num2str(P(:,1)))
%   disp(' ')
%   disp('           peak #    Position      Height     Width        Area')
%   disp(P)
%
% Example 7: Find, measure, and plot noisy peaks with unknown positions
%   x=-50:.2:50;
%   y=exp(-(x).^2)+exp(-(x+50*rand()).^2)+.02.*randn(size(x));
%   plot(x,y,'m.')
%   P=autofindpeaks(x,y);
%   text(P(:,2),P(:,3),num2str(P(:,1)))
%   disp('           peak #    Position    Height       Width        Area')
%   disp(P)
%
% Note: You can pass the detection parameters found by autofindpeaks to 
% other functions, such as measurepeaks. For example:
% x=[0:.1:50];y=5+5.*sin(x)+randn(size(x));[P,A]=autofindpeaks(x,y,3);
% P=measurepeaks(x,y,A(1),A(2),A(3),A(4),1);
%
% Related functions:
% autofindpeaksplot.m, findpeaksG.m, findvalleys.m, findpeaksL.m,
% findpeaksb.m, findpeaksb3.m, findpeaksplot.m, peakstats.m, findpeakSNR.m,
% findpeaksGSS.m, findpeaksSG, findpeaksLSS.m, findpeaksfit.m, findsteps.m,
% findsquarepulse.m, idpeaks.m, measurepeaks.m
% Copyright (c) 2017 Thomas C. O'Haver

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
        disp('Input arguments for findpeaks.... functions: ' )
        disp(['(x,y,'  num2str(SlopeThreshold) ',' num2str(AmpThreshold)  ',' num2str(smoothwidth)  ',' num2str(peakgroup) ',3)'] )
        DetectionParameters=[SlopeThreshold AmpThreshold smoothwidth peakgroup];
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
        disp('Input arguments for findpeaks.... functions: ' )
        disp(['(x,y,'  num2str(SlopeThreshold) ',' num2str(AmpThreshold)  ',' num2str(smoothwidth)  ',' num2str(peakgroup) ',3)'] )
        DetectionParameters=[SlopeThreshold AmpThreshold smoothwidth peakgroup];
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
        disp('Input arguments for findpeaks.... functions: ' )
        disp(['(x,y,'  num2str(SlopeThreshold) ',' num2str(AmpThreshold)  ',' num2str(smoothwidth)  ',' num2str(peakgroup) ',3)'] )
        DetectionParameters=[SlopeThreshold AmpThreshold smoothwidth peakgroup];
end
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
DetectionParameters=[SlopeThreshold AmpThreshold smoothwidth peakgroup];
% Remove comment on next two lines a add plotting
% plot(x,y)
% text(P(:,2),P(:,3),num2str(P(:,1)))
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
%   SegmentedSmooth(y,w,type) smooths with ends=0.
%   SegmentedSmooth(y,w) smooths with type=1 and ends=0.
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

% ----------------------------------------------------------------------
function a=rmnan(a)
% Removes NaNs and Infs from vectors, replacing with nearest real numbers.
% Example:
%  >> v=[1 2 3 4 Inf 6 7 Inf  9];
%  >> rmnan(v)
%  ans =
%     1     2     3     4     4     6     7     7     9
la=length(a);
if isnan(a(1)) || isinf(a(1)),a(1)=0;end
for point=1:la,
    if isnan(a(point)) || isinf(a(point)),
        a(point)=a(point-1);
    end
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