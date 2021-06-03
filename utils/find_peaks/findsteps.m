function P=findsteps(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,peakgroup)
% function P=findsteps(x,y,SlopeThreshold,AmpThreshold,peakgroup)
% Function to locate the positive steps in a noisy x-y time series data
% set. Detects steps by computing the first derivative of y that exceed
% SlopeThreshold, computes the step height by taking the difference between
% the maximum and minimum y values over a number of data point equal to
% "Peakgroup". Returns list (P) containing step number, minimun and maximum
% x and y positions, and the step height of each step detected. Arguments
% "slopeThreshold", "AmpThreshold" and "SmoothWidth" control step
% sensivity. Higher values will neglect smaller features. Increasing
% "SmoothWidth" reduces noise. Related functions:
% findpeaks.m, findvalleys.m, findpeaksL.m, findpeaksb.m, findpeaksplot.m,
% peakstats.m, findpeaksnr.m, findpeaksGSS.m, findpeaksLSS.m,
% findpeaksfit.m,findsquarepulse.m

% Version 2. Copyright (c) 2014, Thomas C. O'Haver
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
peakgroup=round(peakgroup);
y=fastsmooth(y,SmoothWidth,3,1);
d=deriv(y);
P=[0 0 0 0 0 0];
upstep=1;
j=1;
while j<length(y)-peakgroup,
    if d(j) > SlopeThreshold, % if derivative is larger than SlopeThreshold
        startgroup=round(j-peakgroup/2+1);
        endgroup=round(j+peakgroup+1);
        midgroup=round(startgroup+peakgroup/2);
        miny=min(y(startgroup:midgroup));
        maxx=max(x(midgroup:endgroup)); 
        maxy=max(y(midgroup:endgroup));
        basex=x(startgroup);
        basey=y(startgroup);
        % Construct matrix P. One row for each upstep
        if (maxy-miny)>AmpThreshold,
            % count this as a valid peak
            P(upstep,:) = [round(upstep) basex basey maxx maxy maxy-miny];
            upstep=upstep+1; % Move on to next peak
        end %   if (maxy-miny)>AmpThreshold,
%          Uncommment next 3 lines to plot each step, wait for keypress
%          plot(x(startgroup:endgroup),y(startgroup:endgroup))
%          [round(upstep) basex basey maxx maxy maxy-miny]
%          pause
        j=j+peakgroup;
    end
    j=j+1;
end
% ----------------------------------------------------------------------
function d=deriv(a)
% First derivative of vector using 2-point central difference.
%  T. C. O'Haver, 1988.
n=length(a);
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
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
