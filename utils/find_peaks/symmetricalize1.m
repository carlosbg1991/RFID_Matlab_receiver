function ySym=symmetricalize1(t,y,factor,smoothwidth)
% function ySym=symmetricalize(t,y,factor,smoothwidth)
% symmetricalization of exponentially broadened peaks by the weighted
% addition of the first derivative. The area under the peak is not changed.
%
% Example 1: Symmetricalization of an exponentially broadend Gaussian.
% t=[0:.01:20]';
% tau=1; % tau is exponential decay constant in x units 
% y=expgaussian2(t,8,1,tau); 
% factor=1;% For Gaussian, factor is about equal to tau
% symy=symmetricalize(t,y,1,1);
% plot(t,y,t,symy)
%
% Example 2: Symmetricalization of an exponentially broadended Lorentzian.
% t=0:5:3000;
% tau=50; % tau is exponential decay constant in x units 
% y=explorentzian(t,1000,100,-tau);
% factor=250; % For lorentzian, factor is about 5*tau
% symy=symmetricalize(t,y,250,1);
% plot(t,y,t,symy)
% OriginalArea=trapz(t,y)
% SharpenedArea=trapz(t,symy)
%
% Example 3: Symmetricalization of a signal with multiple EMG peaks
% x=DataMatrix3(:,1); % Data file from the SPECTRUM.ZIP archive
% y=DataMatrix3(:,2); 
% tau=17; % tau is exponential decay constant in x units 
% smoothwidth=5;
% symy=symmetricalize(x,y,tau,smoothwidth); 
% plot(x,y,x,symy)

% Copyright (c) 2018 Thomas C. O'Haver
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

FirstDeriv=derivxy(t,fastsmooth(y,smoothwidth,3,1));
% Weighted sum of y and its first derivative
% ySym = symmetricalized y
ySym=y+factor.*FirstDeriv;

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
%
% Examples:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%
% x=1:100;
% y=randn(size(x)); 
% plot(x,y,x,fastsmooth(y,5,3,1),'r')
% xlabel('Blue: white noise.    Red: smoothed white noise.')
%
% Demonstration scripts and functions: SmoothExperiment.m, 
% SmoothWidthTest.m, smoothdemo.m, SmoothOptimization.m, 
% DemoSegmentedSmooth.m, DeltaTest.m, SmoothVsCurvefit.m
%
% Related functions: SegmentedSmooth.m, medianfilter.m, killspikes.m  
% Download from 
% https://terpconnect.umd.edu/~toh/spectrum/functions.html

if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
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
