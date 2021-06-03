function P=findpeaksGw(x,y,SlopeThreshold,AmpThreshold,level,peakgroup)
% function P=findpeaksGw(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Wavelet-based function to locate the positive peaks in a noisy x-y time
% series data set (Requires the "wdenoise" function found in the Wavelet
% toolbox). Detects peaks by looking for downward zero-crossings in the
% first derivative that exceed SlopeThreshold. Returns list (P) containing
% peak number and position, height, width, and area of each peak. Arguments
% "slopeThreshold", "ampThreshold" control peak sensitivity. Level is the
% wavelet level, Higher values will neglect smaller features. "Smoothwidth"
% is the width of the smooth applied before peak detection; larger values
% ignore narrow peaks. If smoothwidth=0, no smoothing is performed.
% "Peakgroup" is the number points around the top part of the peak that are
% taken for measurement. If Peakgroup=0 the local maximum is takes as the
% peak height and position.
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% Version 1. (c) T.C. O'Haver, 2019  
%
% Examples:
% findpeaksGw(0:.01:2,humps(0:.01:2),0,-1,5,5)
% x=[0:.01:50];findpeaksGw(x,cos(x),0,-1,5,5)
% x=[0:.01:5]';findpeaksGw(x,x.*sin(x.^2).^2,0,-1,5,5)
% x=[-10:.1:10];y=exp(-(x).^2);findpeaksGw(x,y,0.005,0.3,3,5);
%
% Find, measure, and plot noisy peak with unknown position
% x=[-10:.2:10];
% y=exp(-(x+5*randn()).^2)+.1.*randn(size(x));
% P=findpeaksGw(x,y,0.003,0.5,6,9);
% xx=linspace(min(x),max(x));
% yy=P(3).*gaussian(xx,P(2),P(4));
% plot(x,y,'.',xx,yy)
%
% Related functions:
% findpeaksG.m, findvalleys.m, findpeaksL.m, findpeaksb.m, findpeaksb3.m,
% findpeaksplot.m, peakstats.m, findpeaksnr.m, findpeaksGSS.m,
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
level=round(level);
peakgroup=round(peakgroup);
if level>1
    d=wdenoise(deriv(y),level);
else
    d=deriv(y);
end
n=round(peakgroup/2+1);
P=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=2*round(level/2)-1:length(y)-level-1
    if sign(d(j)) > sign (d(j+1)) % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold*y(j) % if slope of derivative is larger than SlopeThreshold
            if or(y(j) > AmpTest, y(j+1) > AmpTest)  % if height of peak is larger than AmpThreshold (new version by Anthony Willey)
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup>2
                    [Height, Position, Width]=gaussfit(xx,yy);
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
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold
                    % Skip this peak
                else % Otherwiase count this as a valid peak
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
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1
    d(j)=(a(j+1)-a(j-1)) ./ 2;
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

maxy=max(y);
for p=1:length(y),
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
z=log(y);
coef=polyfit(x,z,2);
a=coef(3);
b=coef(2);
c=coef(1);
Height=exp(a-c*(b/(2*c))^2);
Position=-b/(2*c);
Width=2.35482/(sqrt(2)*sqrt(-c));