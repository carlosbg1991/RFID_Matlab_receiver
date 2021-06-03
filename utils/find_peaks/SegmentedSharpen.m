function SymmetrizedSignal=SegmentedSharpen(y,factor1,factor2,SmoothWidth)
% Segmented peak sharpening (resolution enhancement) function by derivative
% method. The arguments: y is the dependent variable vector;
% factor1 and factor2 are 2nd and 4th derivative weighting factor vectors.
% SmoothWidth is the smooth width vector. Larger values of factor1 and
% factor2 will reduce the peak width but will cause artifacts in the
% baseline near the peak. Adjust the factors for the the best compromise.
% Use minimum smooth width needed to reduce excess noise. Use iSignal to
% adjust these factors. Version 1.2 June 2019. Tom O'Haver toh@umd.edu
% 
% Example: Four equal Gaussians with increasing sharpening applied.
%  x=[0:.01:40];
%  y=exp(-(x-6).^2)+exp(-(x-14).^2)+exp(-(x-22).^2)+exp(-(x-30).^2);
%  factor1=[900 1000 1100 1200];
%  factor2=[1.1e6 1.2e6 1.3e6 1.4e6];
%  smoothwidth=[3 3 3 3];
%  Enhancedsignal=SegmentedSharpen(y,factor1,factor2,smoothwidth);
%  plot(x,y,x,Enhancedsignal,'r')

%  If the peak widths increase regularly across the signal, you can calculate 
%  a reasonable initial value for the factor1, factor2, and SmoothWidth    
%  vectors by giving only the number of segments (NumSegments) and the first 
%  and last values:
%   f1step=(endf1-startf1)/NumSegments;
%   factor1=startf1:f1step:endf1;
%  and similarly for factor2 and SmoothWidth.
%
% Related functions: enhance.m, SegmentedSmooth.m,SegGaussDeconv.m,
% SegExpDeconv.m
 
%   Copyright (c) 2017, Thomas C. O'Haver
%   
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
NumSegments=length(factor1);
SegLength=round(ly./NumSegments);
Enhancedsegment=zeros(ly,NumSegments);
EnhancedSignal=zeros(1,ly);
for Segment=1:NumSegments
    Enhancedsegment(:,Segment)=enhance(y,factor1(Segment),factor2(Segment),SmoothWidth(Segment));
    startindex=(1+(Segment-1)*SegLength);
    endindix=startindex+SegLength-1;
    if endindix>ly,endindix=ly;end
    indexrange=startindex:endindix;
    EnhancedSignal(indexrange)=Enhancedsegment(indexrange,Segment);
end
SymmetrizedSignal=SegmentedSmooth(EnhancedSignal,SmoothWidth,3);

function Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth)
% Resolution enhancement function by derivative method. the
% arguments factor1 and factor 2 are 2nd and 4th derivative weighting
% factors. Larger values of factor1 and factor2 will reduce the 
% peak width but will cause artifacts in the baseline near 
% the peak.  Adjust the factors for the the best compromise. 
% Use minimum smooth width needed to reduce excess noise. 
% Use iSignal to adjust these factors 
% 
% Example:
%  x=[0:.01:18];
%  y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2);
%  Enhancedsignal=enhance(y,1000,1000000,3)
%  plot(x,y,x,Enhancedsignal,'r')
%
d2=secderiv(signal);  % Computes second derivative
d4=secderiv(d2);   % Computes fourth derivative
Enhancedsignal = signal-factor1.*tsmooth(d2,SmoothWidth)+...
factor2.*tsmooth(tsmooth(tsmooth(d4,SmoothWidth),SmoothWidth),SmoothWidth);

function d=secderiv(a)
% Second derivative of vector using 3-point central difference.
% Example: secderiv([0 1 2 3 4 3 2 1 0]) yields [ 0 0 0 0 -2 0 0 0 0]
%  T. C. O'Haver, 2006.
n=length(a);
d=zeros(size(a));
for j = 2:n-1
  d(j)=a(j+1) - 2.*a(j) + a(j-1);
end
d(1)=d(2);
d(n)=d(n-1);

function s=tsmooth(Y,w)
%  tsmooth(Y,w) smooths vector Y by a triangular function of halfwidth w
%  T. C. O'Haver, 1988.
v=ones(1,w);v=conv(v,v);
S=conv(Y,v);
startpoint=(length(v) + 1)/2;
endpoint=length(Y)+startpoint-1;

s=S(startpoint:endpoint) ./ sum(v);

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

if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
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