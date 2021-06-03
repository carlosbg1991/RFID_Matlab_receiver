function SmoothSymmSignal=symmetricalize(x,y,factor,smoothwidths,type,ends)
% function ySym=symmetricalize(t,y,factor,smoothwidth,type,ends)
% Segmented constant-area symmetricalization of exponentially broadened
% peaks by the weighted addition of the first derivative. Divides y into a
% number of equal-length segments defined by the length of the vector
% 'factor', then adds each segment to its first derivative times the value
% of "factor" for that segment. The "strength" of the effect is determined
% by "factor", which can be a scalar or a vector. The correct value(s) of
% "factor" depends on the peak shape; for an exponentially modified
% Gaussian it's equal to the exponential time constant tau (1/lambda in
% some formulations of the EMG), and for an exponentially broadened
% Lorentzian, it's about 5 times tau. If the 4th input argument is
% supplied, the result is then smoothed by the segmented smooth widths
% determined by "smoothwidths", which can also be a scalar or vector (not
% necessarily of the same length as "factor"). Version 1.1, 2018.
%
% Example 1: Symmetricalization of an exponentially broadend Gaussian.
% t=[0:.01:20]';
% tau=1; % tau is exponential decay constant in x units (1 segment)
% y=expgaussian2(t,8,1,tau); 
% factor=1;% For Gaussian, factor is about equal to tau
% symy=symmetricalize(t,y,1);
% plot(t,y,t,symy)
% xlabel('x')
% ylabel('y')
% title('Blue: Original signal    Red: Symmetrized by first derivative addition') 
%
% Example 2: Symmetricalization of an exponentially broadended Lorentzian.
% t=0:5:3000;
% tau=50; % tau is exponential decay constant in x units (1 segment)
% y=explorentzian(t,1000,100,-tau);
% factor=250; % For lorentzian, factor is about 5*tau
% symy=symmetricalize(t,y,250);
% plot(t,y,t,symy)
% xlabel('x')
% ylabel('y')
% title('Blue: Original signal    Red: Symmetrized by first derivative addition') 
% OriginalArea=trapz(t,y)
% SharpenedArea=trapz(t,symy)
% xlabel('x')
% ylabel('y')
% title('Blue: Original signal    Red: Symmetrized by first derivative addition') 
%
% Example 3: Symmetricalization of a signal with multiple EMG peaks
% load DataMatrix3; % Data file from the SPECTRUM.ZIP archive
% x=DataMatrix3(:,1);
% y=DataMatrix3(:,2); 
% tau=17; % tau is exponential decay constant in x units (1 segment)
% smoothwidth=5;
% symy=symmetricalize(x,y,tau,smoothwidth); 
% plot(x,y,x,symy)
% xlabel('x')
% ylabel('y')
% title('Blue: Original signal    Red: Symmetrized by first derivative addition') 
%
% Example 4: Segmented symmetricalization of 5 Gaussians with different
% exponential broadening time constants.
% x=5:.1:65;
% y=modelpeaks2(x, [1 5 5 5 5], [1 1 1 1 1], [10 20 30 40 50], [3 3 3 3 3], [0 -5 -10 -15 -20]);
% factor=[.4 .5 .6 1.1 1.3 1.5 2 2]; % Defines 8-segment symmetricalization factor
% Symmy=symmetricalize(x,y,factor,[3 5],1,1);
% plot(x,y,x,Symmy)
% xlabel('x')
% ylabel('y')
% title('Blue: Original signal    Red: Symmetrized by first derivative addition') 
%
% Copyright (c) 2019 Thomas C. O'Haver
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
if nargin==3, smoothwidths=1; ends=0; type=1; end
if nargin==4, ends=0; type=1; end
if nargin==5, ends=0; end
ly=length(y);
NumSegments=length(factor);
SegLength=round(ly./NumSegments);
SymmSegment=zeros(ly,NumSegments);
SymmSignal=zeros(1,ly);
for Segment=1:NumSegments
    SymmSegment(:,Segment)=y+factor(Segment).*derivxy(x,y);
    startindex=(1+(Segment-1)*SegLength);
    endindix=startindex+SegLength-1;
    if endindix>ly,endindix=ly;end
    indexrange=startindex:endindix;
    SymmSignal(indexrange)=SymmSegment(indexrange,Segment);
end
SmoothSymmSignal=SegmentedSmooth(SymmSignal,smoothwidths,type,ends);

function d=derivxy(x,y)
% First derivative of y with respect to x using 2-point central difference.
%  T. C. O'Haver, 2011. Corrected 2/4/14
% Example:
% x = [1 2 4 7 10 14];y = 2*x;derivxy(x,y)
% ans =
%      2     2     2     2     2     2
n=length(y);
d=zeros(size(y));
d(1)=(y(2)-y(1))./(x(2)-x(1));
d(n)=(y(n)-y(n-1))./(x(n)-x(n-1));
for j = 2:n-1
  d(j)=(y(j+1)-y(j-1)) ./ ((x(j+1)-x(j-1)));
end

function SmoothedSignal=SegmentedSmooth(y,smoothwidths,type,ends)
%   SegmentedSmooth(y,w,type,ends) divides y into a number of equal-length
%   segments defined by the length of the vector 'smoothwidths', then
%   smooths each segment with smooth of type 'type' and width defined by
%   the elements of vector 'smoothwidths'. Version 1, November 2016.
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

%   Example 1: Three-segment smooth of random white noise, 
%   smooth widths of 2,20, and 200:
%      x=1:10000;y=randn(size(x));
%      plot(x,SegmentedSmooth(y,[2 20 200],3,0))
%
%   Example 2: 20-segment smooth, odd smooth widths from 1 to 41,  
%   calculated by colon notation:
%      plot(x,SegmentedSmooth(y,[1:2:41],3,0))
%
% Related functions: fastsmooth.m, SegmentedSharpen.m,SegGaussDeconv.m,
% SegExpDeconv.m

%   Copyright (c) 2016, Thomas C. O'Haver
   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in
%   all copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
%   THE SOFTWARE.

ly=length(y);
NumSegments=length(smoothwidths);
SegLength=round(ly./NumSegments);
SmoothSegment=zeros(ly,NumSegments);
SmoothedSignal=zeros(1,ly);
for Segment=1:NumSegments
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
for k=1:L-w
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
  if ends==1
  startpoint=(smoothwidth + 1)/2;
  SmoothY(1)=(Y(1)+Y(2))./2;
  for k=2:startpoint
     SmoothY(k)=mean(Y(1:(2*k-1)));
     SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
  end
  SmoothY(L)=(Y(L)+Y(L-1))./2;
  end