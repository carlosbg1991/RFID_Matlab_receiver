function P=findpeaksG2d(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% P=findpeaksG2d(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype) 
% Function to locate the positive peaks and shoulders in a noisy x-y time
% series data set. Detects peaks in the negative of the second derivative
% of the signal, by looking for downward slopes in the third derivative
% that exceed SlopeThreshold. Returns list (P) containing peak number and
% position and an estimate of the height, width, and area of each peak. For
% more accurate measurement, use curve fitting (peakfit.m).
%
% Input arguments "slopeThreshold", "ampThreshold" and "smoothwidth"
% control peak sensitivity. Higher values will neglect smaller features.
% "Smoothwidth" is the width of the smooth applied before peak detection;
% larger values ignore narrow peaks. If smoothwidth=0, no smoothing is
% performed. "Peakgroup" is the number points around the top part of the
% peak that are taken for measurement. If Peakgroup=0 the local maximum is
% takes as the peak height and position. 
% "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar)
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
%
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver, December 2016. 
%
% Example: run TestFindpeaksG2d.m, which shows that this function can
% locate peaks resultin in 'shoulders' that  do not produce a distinct
% maximum in the original signal.
%
% Related functions:
% findpeaksSG.m, findvalleys.m, findpeaksL.m, findpeaksb.m, findpeaksb3.m,
% findpeaksplot.m, peakstats.m, findpeaksnr.m, findpeaksGSS.m, ipeak.m
% findpeaksLSS.m, findpeaksfit.m, findsteps.m, findsquarepulse.m, idpeaks.m

% Copyright (c) 2013, 2014 Thomas C. O'Haver
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
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=-fastsmooth(deriv3(y),smoothwidth,smoothtype);
%    figure(2);plot(x,d)
else
    d=y;
end
n=round(peakgroup/2+1);
P=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth-1,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold, % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpThreshold,  % if height of peak is larger than AmpThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));dd=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    % [min(groupindex) max(groupindex)]
                    xx(k)=x(groupindex);
                    yy(k)=y(groupindex);
                    dd(k)=d(groupindex);
                end % for k=1:peakgroup,...
                if peakgroup>2,
                    [Height, Position, Width]=gaussfit(xx,dd);
                    PeakX=real(Position);   % Compute peak position and height of fitted parabola
                    PeakY=max(yy);
                    MeasuredWidth=real(Width);
                    % if the peak is too narrow for least-squares technique to work
                    % well, just use the max value of y in the sub-group of points near peak.
                else
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end % if peakgroup>2,...
                % Construct matrix P. One row for each peak detected,
                % containing the peak number, peak position (x-value) and
                % peak height (y-value). If peak measurement fails and
                % results in NaN, or if the measured peak height is less
                % than AmpThreshold, skip this peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold,
                    % Skip this peak
                else % Otherwise count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth];
                    peak=peak+1; % Move on to next peak 
                end % if isnan...
            end % if y(j) > AmpThreshold,...
        end %  if d(j)-d(j+1) >...
    end %  if sign(d(j)) >...
end % for j=2*round(smoothwidth/2)-1:...
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

function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth
%  of width w. Version 2.0, May 2008.
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar)
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively
%     smaller smooths the closer to the end. (In this mode the
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%  T. C. O'Haver, May, 2008.
if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
switch type
    case 1
        SmoothY=sa(Y,w,ends);
    case 2
        SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
        SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
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
% [Height,Position,Width]=gaussfit(x,y)
% plot(x,y,'o',linspace(0,8),Height.*gaussian(linspace(0,8),Position,Width))
% Copyright (c) 2012, Thomas C. O'Haver

% To prevent problems from taking the log of zero or negative values,
% make the lowest value of y equal to 1% of the maximum value.
maxy=max(y);
for p=1:length(y),
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % end of for p=1:length(y),
    logyyy=log(abs(y));
    [coef,S,MU]=polyfit(x,logyyy,2);
    c1=coef(3);c2=coef(2);c3=coef(1);
    % Compute peak position and height or fitted parabola
    Position=-((MU(2).*c2/(2*c3))-MU(1));
    Height=exp(c1-c3*(c2/(2*c3))^2);
    Width=norm(MU(2).*2.35703/(sqrt(2)*sqrt(-1*c3)));
    
function d=deriv2(a)
% Second derivative of vector using 3-point central difference.
%  T. C. O'Haver, 2006.
n=length(a);
for j = 2:n-1;
  d(j)=a(j+1) - 2.*a(j) + a(j-1);
end
d(1)=d(2);
d(n)=d(n-1);