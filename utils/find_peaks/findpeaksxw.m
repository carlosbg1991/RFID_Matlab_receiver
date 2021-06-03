function P=findpeaksxw(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% function P=findpeaksxw(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Function to locate the positive peaks in a noisy x-y time series data
% set. Detects peaks by looking for downward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (P) containing peak number, position, height, and width
% (full width at half maximum) of each peak. Arguments "slopeThreshold",
% "ampThreshold" and "smoothwidth" control peak sensitivity. Higher values
% will neglect smaller features. "Smoothwidth" is the width of the smooth
% applied before peak detection; larger values ignore narrow peaks. If
% smoothwidth=0, no smoothing is performed. "Peakgroup" is the number
% points around the top part of the peak that are taken for measurement. If
% the value of PeakGroup is 1 or 2, the maximum y value of the 1 or 2
% points at the point of zero-crossing is taken as the peak height value;
% if PeakGroup is n is 3 or greater, the average of the next n points is
% taken as the peak height value. For spikes or very narrow peaks, keep
% PeakGroup=1 or 2; for broad or noisy peaks, make PeakGroup larger to
% reduce the effect of noise. The argument "smoothtype" determines the
% smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar)
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver.  Version 1, March, 2018  
%
% Examples:
% x=0:.01:2;findpeaksxw(x,humps(x).^2,0,-1,1,1)
%
% Regularly spaced peaks of constant height and width
% x=[0:.01:50];findpeaksxw(x,cos(x),0,-1,5,5)
%
% Peaks get taller but narrower as x increases:
% x=[0:.01:5]';y=x.*sin(x.^2).^2;plot(x,y);findpeaksxw(x,y,0,-1,5,5)
%
% Find, plot, and number noisy peaks with unknown positions
% x=-50:.2:50;
% y=exp(-(x).^2)+exp(-(x+50*rand()).^2)+exp(-(x+50*rand()).^2)+.02.*randn(size(x));
% plot(x,y,'m.')
% P=findpeaksxw(x,y,0.001,0.2,1,2,3);
% text(P(:,2),P(:,3),num2str(P(:,1)))
% disp('    peak #    Position   Height   FWHM')
% disp(P)
%
% Related functions:
% findpeaksx.m, findvalleys.m, findpeaksG.m, findpeaksb.m, findpeaksb3.m,
% findpeaksplot.m, peakstats.m, findpeaksnr.m, findpeaksGSS.m,
% findpeaksLSS.m, findpeaksfit.m, findsteps.m, findsquarepulse.m, idpeaks.m

% Copyright (c) 2013, 2018 Thomas C. O'Haver
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
if smoothwidth>1
    d=fastsmooth(deriv(y),smoothwidth,smoothtype);
else
    d=deriv(y);
end
n=round(peakgroup/2+1);
P=[0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
PeakX=0;
PeakY=0;
% plot(x,d,'.')
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth-1
    if sign(d(j)) > sign (d(j+1)) % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold % if slope of derivative is larger than SlopeThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup<3
                    PeakY=max(yy);
                else
                    PeakY=mean(yy);
                end
                pindex=val2ind(yy,PeakY);
                PeakX=xx(pindex(1));
            % Construct matrix P. One row for each peak detected,
            % containing the peak number, peak position (x-value),
            % peak height (y-value), and halfwidth. If peak measurement fails and
            % results in NaN, or if the measured peak height is less
            % than AmpThreshold, skip this peak
            if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold
                % Skip this peak if PeakY is below AmpThreshold
            else % Otherwise count this as a valid peak
                P(peak,:) = [round(peak) PeakX PeakY halfwidth(x,y,PeakX)];
                peak=peak+1; % Move on to next peak
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

function [FWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,y,xo)
% function [FWHM,slope1,slope2,x1,x2]=halfwidth(x,y,xo) computes the full
% width at half maximum of the maximum of any shape peak that peaks near xo
% and has a zero baseline. If xo is omitted, it computes the halfwidth from
% the maximum y value. Not highly accurate if the function is too sparsely
% sampled. Also optionally returns the leading and trailing edge slopes
% (slope1 and slope 2, respectively) and the leading and trailing edge
% half-width-at-half-maximum, if those output arguments are incuded. 
% Tom O'Haver (toh@umd.edu) Version 2, July 2016
%
% Example 1:
% x=-5:.01:5;
% y=sinc(x).^2;
% FWHM=halfwidth(x,y,0) % Area of central peak
% FWHM=halfwidth(x,y,1.5) % Area of positive side peak
%
% Example 2:
% x=[0:.1:10];
% W=3; % W is the true half-width
% y=gaussian(x,5,W);
% FWHM=halfwidth(x,y,5)
% FWHM=halfwidth(x,y) % gives the same result since there is only one peak
%
% Example 3: % Returns NaN because signal does not contain the halfwidth
% x=[0:.1:1];
% y=exp(-(x).^2)
% halfwidth(x,y) 
% 
% Example 4: Also returns the leading and trailing edge slopes (+1 and -1)
% and the leading and trailing edge half-width-at-half-maximum
% x=[0:.1:10];
% y=triangle(x,5,1);
% [FWHM,slope1,slope2,x1,x2]=halfwidth(x,y)
%
% Example 5. Comparison of Gaussian (g) and BiGaussian peaks (bg)
% x=[1:100];
% g=gaussian(x,50,20);
% bg=bigaussian(x,50,20,3);
% plot(x,g,x,bg)
% [GaussianFWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,g)
% [BiGaussianFWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,bg)

if nargin==2,xo=x(val2ind(y,max(y)));end
try   
    indmax=val2ind(x,xo);
    maxy=y(indmax);
    oy=y-maxy/2;
    
    n=indmax;
    while oy(n)>0
        n=n-1;
    end
    x1=interp1([y(n) y(n+1)],[x(n) x(n+1)],maxy/2);
    slope1=(y(n+1)-y(n))./(x(n+1)-x(n));
    
    n=indmax;
    while oy(n)>0
        n=n+1;
    end
    x2= interp1([y(n-1) y(n)],[x(n-1) x(n)],maxy/2);
    slope2=(y(n-1)-y(n))./(x(n-1)-x(n));
    
    FWHM=x2-x1;
    hwhm1=xo-x1;
    hwhm2=x2-xo;
catch
    FWHM=NaN;
end
