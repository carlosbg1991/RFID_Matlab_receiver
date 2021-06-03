function EMLplusfirstderivative
% Demonstration of symmetricalization of an exponentially modified
% Lorentzian (EML) by weighted addition of the first derivative, weighting
% factor=tau=1/lambda. Shows that the peak width is reduced, the
% peak made more symmetrical, and the area is preserved.
s=100;
lambda=.02;
mu=1000;
deltat=1;
Noise=0;
t=0:deltat:2500;
figure(1);
clf
tau=1/lambda;

factor=tau;
EML=explorentzian(t,mu,s,-tau)+Noise.*randn(size(t))';
FirstDeriv=derivxy(t,EML);
% Weighted sum of EML and its first derivative
% EMLSym = symmetricalized EML
EMLSym=EML+factor.*FirstDeriv;
% Fit to Lorentzian to show how close the result is to a Gaussian
[~,GOF]= peakfit([t;EMLSym'],0,0,1,2,0,1,0,0,0,0);
Rsquared=GOF(2);
[HalfwidthBefore,slopebefore1,slopebefore2]=halfwidth(t,EML,mu);
SlopeRatioBefore=slopebefore1/slopebefore2;
[HalfwidthAfter,slopeafter1,slopeafter2]=halfwidth(t,EMLSym,mu);
SlopeRatioAfter=slopeafter1/slopeafter2;
disp(' ')
disp(['Halfwidth of original Lorentzian: ' num2str(s)])
disp(['Halfwidth of EML: ' num2str(HalfwidthBefore)])
disp(['Slope Ratio Before: ' num2str(SlopeRatioBefore)])
disp(['Area Before: ' num2str(trapz(t,EML))])
disp(' ')
disp(['First derivative weight: ' num2str(factor)])
disp(['Halfwidth After: ' num2str(HalfwidthAfter)])
disp(['Slope Ratio After: ' num2str(SlopeRatioAfter)])
disp(['R2 for fit to Lorentzian: ' num2str(Rsquared)])
disp(['Area After: ' num2str(trapz(t,EMLSym))])
plot(t,EML,t,EMLSym,t,factor.*derivxy(t,EML),'.');
title('First derivative symmetricalization of an exponentially modified Lorentzian (EML)')
xlabel('Blue: Original EML.   Yellow: 1st derivative * weight.    Red: EML plus 1st derivative * weight')
text(0,.5,['  Width (sigma) = ' num2str(s)])
text(0,.4,['  Asymmetry (lambda) = ' num2str(lambda)])
text(0,.3,['  Asymmetry (tau) = ' num2str(tau)])
text(0,.2,['  Weighting factor = ' num2str(factor)])

disp(' ')
% % Additional 2nd and 4th derivative sharpening of the symmetricalized peak, EMLSym
% figure(2)
%  clf 
%  k1=((HalfwidthAfter/deltat).^2)/5;
%  k2=((HalfwidthAfter/deltat).^4)/500;
%  SmoothWidth=1;
%  sharpy=enhance(EMLSym,k1,k2,SmoothWidth);
%  plotrange=1:length(t);
% plot(t(plotrange),EMLSym(plotrange),t(plotrange),sharpy(plotrange))
% title('Further 2nd and 4th derivative sharpening of the symmetricalized peak')
% xlabel('t')

% disp(['Halfwidth after additional sharpening: ' num2str(halfwidth(t,sharpy,mu))])
 
function g = explorentzian(x,pos,wid,timeconstant)
%  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2013
g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
g = ExpBroaden(g',timeconstant);

function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant "t" by multiplying Fourier transforms and inverse
% transforming the result.
% Example:
% clg;
% x=1:100
% y=zeros(size(x));
% y(40:60)=1;
% plot(x,y,x,ExpBroaden(y',-5))
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);

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