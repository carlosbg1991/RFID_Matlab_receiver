function [EMGSym,tau]=AutoSymmetrize(t,EMG,SmoothWidth,plots)
% [EMGSym,tau]=AutoSymmetrize(t,EMG,SmoothWidth,plots)
% Function that automatically determines the exponential time constant tau
% of an exponentially broadened peak defined by [t,EMG], using the first
% derivative addition method. "SmoothWidth" provides for smoothing noise
% signals. If plots=0, not graphs are plotted. If plots=1, three plots
% describing the method are propared. Returns the symmetrized peak EMGSym
% and two estimates of the tau, (1) by zero baseline and (2) by symmetry. 
% Tom O'Haver, July 2019
% Example: 
% A=1;
% s=100; % standard deviation of Gaussian peak, in time units
% tau=100; % Exponential broadening factor, in time units
% mu=1000; % peak center, in time units
% deltat=10; % sampling interval, in time units
% Noise=.001; % random white noise
% plots=1;  % 1 = show plots, 0 = no plots
% SmoothWidth=5; % Smooth width to reduce effect of noise, in points
% t=0:deltat:2500; % time axis
% EMG=A.*expgaussian2(t,mu,s,tau)+Noise.*randn(size(t));
% [EMGSym,tau]=AutoSymmetrize(t,EMG,SmoothWidth,plots);

EMG=fastsmooth(EMG,SmoothWidth,2);
mu=t(val2ind(EMG,max(EMG)));
[hw,slopeafter1,slopeafter2]=halfwidth(t,EMG,mu);
r=slopeafter1/slopeafter2;
EstimatedTau=r.*-hw/4-hw/8;
FirstDeriv=derivxy(t,EMG);
trial=0;
StartFactor=EstimatedTau/2;
EndFactor=EstimatedTau*2;
for factor=StartFactor:1:EndFactor
   trial=trial+1;
   % Weighted sum of EMG and its first derivative
   % EMGSym = symmetricalized EMG
   EMGSym=EMG+factor.*FirstDeriv;
   EMGSym=fastsmooth(EMGSym,SmoothWidth,2);
   Miny(trial)=min(EMGSym);
   TrialTau(trial)=factor;
   MinD=min(FirstDeriv);
   EMGmin(trial)=EMGSym(val2ind(FirstDeriv,MinD));
   [~,slopeafter1,slopeafter2]=halfwidth(t,EMGSym,mu);
   SlopeRatio(trial)=slopeafter1/slopeafter2;
end
polyorder=3;
coef=polyfit(SlopeRatio,TrialTau,polyorder);
tau=polyval(coef,-1); % TrialTau=tau when SlopeRatio equals -1
EMGSym=EMG+tau.*FirstDeriv; % Symmetrize with optimum tau
EMGSym=fastsmooth(EMGSym,SmoothWidth,2);

if plots
    figure(2)
    clf
    plotit(SlopeRatio,TrialTau,polyorder);
    grid
    ylabel('TrialTau (time units)')
    xlabel('Ratio of leading edge to trailing edge slope')
    title(['Ratio equals -1 when peak is symmetrical, at tau = '  num2str(tau) '.'])

    figure(1);
    clf

    plot(t,EMG,t,EMGSym)
    xlabel('Blue: Original exponentially broadened peak.   Red: Symmetrized peak')
    title(['Automated determination of exponential broadening factor, tau = '  num2str(tau) '.'])
    disp(' ')
end

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

function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function s=tsmooth(Y,w)
%  tsmooth(Y,w) smooths vector Y by a triangular function of halfwidth w
%  T. C. O'Haver, 1988.
v=ones(1,w);v=conv(v,v);
S=conv(Y,v);
startpoint=(length(v) + 1)/2;
endpoint=length(Y)+startpoint-1;

s=S(startpoint:endpoint) ./ sum(v);

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
%
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

% Copyright (c) 2012,2018 Thomas C. O'Haver
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


% Copyright (c) 2018, Thomas C. O'Haver
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