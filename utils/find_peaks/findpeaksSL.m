function P=findpeaksSL(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Segmented peak finder, same syntax as findpeaksL except the 3rd to 6th
% input arguments can be vectors with one entry for each segment.
% function P=findpeaksSL(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
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
%   y=lorentzian(x,20,1.5)+lorentzian(x,80,30)+.02.*randn(size(x));
%   plot(x,y,'c.')
%   P=findpeaksSL(x,y,[0.001 .0001],[.2 .2],[5 10],[10 50],3)
%   text(P(:,2),P(:,3),num2str(P(:,1)))
%   disp(' ')
%   disp('   peak #    Position      Height     Width      Area')
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
                    z=rmnan(ones(size(xx))./yy);
                    coef=polyfit(xx,z,2);
                    PeakY=4*coef(1)./((4*coef(1)*coef(3))-coef(2)^2);
                    PeakX=-coef(2)/(2*coef(1));
                    MeasuredWidth=sqrt(abs(((4*coef(1)*coef(3))-coef(2)^2)./coef(1))./sqrt(abs(coef(1))));                    % if the peak is too narrow for least-squares technique to work
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


