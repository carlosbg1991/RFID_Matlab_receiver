function rG=DEMSymm(t,DEMG,L1guess,L2guess,DeltaGuess)
% function
% rG=DEMSymm(t,DEMG,L1guess,L2guess,DeltaGuess)
% Function for the symmetrization of a double exponentially modified peaks
% (DEMP) by two successive applications of weighted addition of the first
% derivative, with weighting factors ideally equal to the two taus. The
% objective is to make the peaks more symmetrical and narrower while
% preserving the peak area. Uses a three-level plus-and-minus 1% bracketing
% technique to help you to determine the best values for the
% first-derivative weighting factors L1 and L2. See the demo script
% Demo1stDerivManualSymm.m. 
%
% Example. Here is a series of command line calls for a typical manual
% search for the best weighting factors, starting with very poor initial
% values. Start with L1=0 and adjust L2, then adjust L1. Figure(2) shows an
% expanded y scale so you can see small ajdustments.
% Load a data file containing x,y, values for a double exponential
% Gaussian.
% >> load DEMG
% % Start with L1=0 and adjust L2:
% >> rG=DEMSymm(t,DEMG,0,200,1);
% >> rG=DEMSymm(t,DEMG,0,300,1);
% >> rG=DEMSymm(t,DEMG,0,400,1);
% >> rG=DEMSymm(t,DEMG,0,350,1);
% >> rG=DEMSymm(t,DEMG,0,340,1);
% >> rG=DEMSymm(t,DEMG,0,330,1);
% >> rG=DEMSymm(t,DEMG,0,335,1);
% >> rG=DEMSymm(t,DEMG,0,333,1);
% % Now vary L1:
% >> rG=DEMSymm(t,DEMG,50,333,1);
% >> rG=DEMSymm(t,DEMG,70,333,1);
% >> rG=DEMSymm(t,DEMG,80,333,1);
% >> rG=DEMSymm(t,DEMG,90,333,1);
% >> rG=DEMSymm(t,DEMG,100,333,1);
% >> rG=DEMSymm(t,DEMG,102,333,1);
% >> rG=DEMSymm(t,DEMG,101,333,1);
% >> rG=DEMSymm(t,DEMG,100,333,1);
% >> rG=DEMSymm(t,DEMG,100.5,333,1);
% The true values of tau1 and tau2 were 100 and 333.3 in this case, so
% these results are quite close.
figure(1)
clf
DEMGPFD=L2guess.*derivxy(t',DEMG);
plot(t,DEMG,'b',t,DEMG+DEMGPFD,'g');
hold on
DEMGPFD=(L2guess+L2guess*DeltaGuess/10).*derivxy(t',DEMG);
plot(t,DEMG+DEMGPFD,'g');
DEMGPFD=(L2guess-L2guess*DeltaGuess/10).*derivxy(t',DEMG);
plot(t,DEMG+DEMGPFD,'g');
DEMGPFD=L2guess.*derivxy(t',DEMG);

CorrDEMG=DEMG+DEMGPFD;
rG=CorrDEMG+L1guess.*derivxy(t',CorrDEMG);
plot(t,DEMG,t,DEMG+DEMGPFD,'g',t,rG,'r');
hold on
CorrDEMG=DEMG+DEMGPFD;
rG=CorrDEMG+(L1guess+L1guess*DeltaGuess/10).*derivxy(t',CorrDEMG);
plot(t,DEMG,'b',t,DEMG+DEMGPFD,'g',t,rG,'r');
CorrDEMG=DEMG+DEMGPFD;
rG=CorrDEMG+(L1guess-L1guess*DeltaGuess/10).*derivxy(t',CorrDEMG);
plot(t,DEMG,'b',t,DEMG+DEMGPFD,'g',t,rG,'r');
grid
hold off
rG=CorrDEMG+L1guess.*derivxy(t',CorrDEMG);
ylabel('y')
xlabel('x')
title('Blue: Original double exponential peak   Green: second stage   Red: Recovered peak')

startx=min(t);
endx=max(t);
miny=min(rG);
maxy=max(rG);
figure(2)
clf
DEMGPFD=L2guess.*derivxy(t',DEMG);
plot(t,DEMG+DEMGPFD,'g');axis([startx endx -maxy/200 maxy/70])
hold on
DEMGPFD=(L2guess+L2guess*DeltaGuess/10).*derivxy(t',DEMG);
plot(t,DEMG+DEMGPFD,'g');axis([startx endx -maxy/200 maxy/70])
DEMGPFD=(L2guess-L2guess*DeltaGuess/10).*derivxy(t',DEMG);
plot(t,DEMG+DEMGPFD,'g');axis([startx endx -maxy/200 maxy/70])
DEMGPFD=L2guess.*derivxy(t',DEMG);

CorrDEMG=DEMG+DEMGPFD;
rG=CorrDEMG+L1guess.*derivxy(t',CorrDEMG);
plot(t,DEMG,'b',t,DEMG+DEMGPFD,'g',t,rG,'r');axis([startx endx -maxy/200 maxy/70])
hold on
CorrDEMG=DEMG+DEMGPFD;
rG=CorrDEMG+(L1guess+L1guess*DeltaGuess/10).*derivxy(t',CorrDEMG);
plot(t,DEMG,'b',t,DEMG+DEMGPFD,'g',t,rG,'r');axis([startx endx -maxy/200 maxy/70])
CorrDEMG=DEMG+DEMGPFD;
rG=CorrDEMG+(L1guess-L1guess*DeltaGuess/10).*derivxy(t',CorrDEMG);
plot(t,DEMG,'b',t,DEMG+DEMGPFD,'g',t,rG,'r');axis([startx endx -maxy/200 maxy/70])
grid
title('Blue: Original double exponential peak   Green: second stage   Red: Recovered peak')
hold off
ylabel('Expanded scale')
% Return rG with middle value of L1
rG=CorrDEMG+L1guess.*derivxy(t',CorrDEMG);


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
for j = 2:n-1;
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
% Related functions: SegmentedSharpen.m

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