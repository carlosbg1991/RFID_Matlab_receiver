% function EMGplusfirstderivative
% Demonstration of symmetricalization of an exponentially modified
% Gaussian (EMG) by weighted addition of the first derivative, weighting
% factor equals tau (the exponential time cpnstant). Shows that the peak
% width is reduced, the peak made more symmetrical, and the area is
% preserved. Automatically determines the optimum tau two different ways,
% by (1) zero baseline and (2) by symmetry. The estimate by symmetry is
% more robust in the presence of noise.
format compact
format short g
clear
A=1;
s=100; % standard deviation of Gaussian peak
tau=200; % Exponential broadening factor
mu=1000; % peak center
deltat=5; % sampling interval
Noise=0.01; % random white noise
SmoothWidth=51; % Smooth width to reduce effect of noise
t=0:deltat:2500; % time axis
EMG=A.*expgaussian2(t,mu,s,tau)+Noise.*randn(size(t));
EMG=fastsmooth(EMG,SmoothWidth,2);

[hw,slopeafter1,slopeafter2]=halfwidth(t,EMG,mu);
r=slopeafter1/slopeafter2;
EstimatedTau=r.*-hw/4-hw/8;
mu=t(val2ind(EMG,max(EMG)));
FirstDeriv=derivxy(t,EMG);
% InitialWidth=halfwidth(t,EMG,mu);
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
   SlopeRatioAfter(trial)=slopeafter1/slopeafter2;
   % disp([trial factor EMGmin(trial) SlopeRatioAfter(trial)])
end

coef=polyfit(SlopeRatioAfter,TrialTau,3);
EstimateBySymmetry=polyval(coef,-1);

for tt=length(TrialTau):-1:1
    if Miny(tt)>.1*min(Miny)
        EstimateByZeroBaseline=TrialTau(tt);
        break
    end
end

figure(2)
clf
coef=plotit(SlopeRatioAfter,TrialTau,3);
grid
ylabel('TrialTau (time units)')
xlabel('Ratio of leading edge to trailing edge slope')
title(['Ratio equals -1 when peak is symmetrical, at tau = '  num2str(EstimateBySymmetry) '.'])

figure(3)
clf
plot(TrialTau,Miny)
xlabel('TrialTau (time units)')
 title(['Graph hits zero at tau = '  num2str(EstimateByZeroBaseline) '.'])
ylabel('Minimum of EMGSym')


figure(1);
clf
factor=EstimateBySymmetry;
EMGSym=EMG+factor.*FirstDeriv;
EMGSym=fastsmooth(EMGSym,SmoothWidth,2);
plot(t,EMG,t,EMGSym)
xlabel('Blue: Original exponentially broadened peak.   Red: Symmetrized peak')
title('Automated determination of exponential broadening factor, tau')
disp('        Noise      SmoothWidth      tau    Symmetry    ZeroBaseline ')
disp([Noise SmoothWidth tau EstimateBySymmetry EstimateByZeroBaseline])
% 
% function d=derivxy(x,y)
% % First derivative of y with respect to x using 2-point central difference.
% %  T. C. O'Haver, 2011. Corrected 2/4/14
% % Example:
% % x = [1 2 4 7 10 14];y = 2*x;derivxy(x,y)
% % ans =
% %      2     2     2     2     2     2
% n=length(y);
% d=zeros(size(y));
% d(1)=(y(2)-y(1))./(x(2)-x(1));
% d(n)=(y(n)-y(n-1))./(x(n)-x(n-1));
% for j = 2:n-1;
%   d(j)=(y(j+1)-y(j-1)) ./ ((x(j+1)-x(j-1)));
% end
% 
% function [FWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,y,xo)
% % function [FWHM,slope1,slope2,x1,x2]=halfwidth(x,y,xo) computes the full
% % width at half maximum of the maximum of any shape peak that peaks near xo
% % and has a zero baseline. If xo is omitted, it computes the halfwidth from
% % the maximum y value. Not highly accurate if the function is too sparsely
% % sampled. Also optionally returns the leading and trailing edge slopes
% % (slope1 and slope 2, respectively) and the leading and trailing edge
% % half-width-at-half-maximum, if those output arguments are incuded. 
% % Tom O'Haver (toh@umd.edu) Version 2, July 2016
% %
% % Example 1:
% % x=-5:.01:5;
% % y=sinc(x).^2;
% % FWHM=halfwidth(x,y,0) % Area of central peak
% % FWHM=halfwidth(x,y,1.5) % Area of positive side peak
% %
% % Example 2:
% % x=[0:.1:10];
% % W=3; % W is the true half-width
% % y=gaussian(x,5,W);
% % FWHM=halfwidth(x,y,5)
% % FWHM=halfwidth(x,y) % gives the same result since there is only one peak
% %
% % Example 3: % Returns NaN because signal does not contain the halfwidth
% % x=[0:.1:1];
% % y=exp(-(x).^2)
% % halfwidth(x,y) 
% % 
% % Example 4: Also returns the leading and trailing edge slopes (+1 and -1)
% % and the leading and trailing edge half-width-at-half-maximum
% % x=[0:.1:10];
% % y=triangle(x,5,1);
% % [FWHM,slope1,slope2,x1,x2]=halfwidth(x,y)
% %
% % Example 5. Comparison of Gaussian (g) and BiGaussian peaks (bg)
% % x=[1:100];
% % g=gaussian(x,50,20);
% % bg=bigaussian(x,50,20,3);
% % plot(x,g,x,bg)
% % [GaussianFWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,g)
% % [BiGaussianFWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,bg)
% 
% if nargin==2,xo=x(val2ind(y,max(y)));end
% try   
%     indmax=val2ind(x,xo);
%     maxy=y(indmax);
%     oy=y-maxy/2;
%     
%     n=indmax;
%     while oy(n)>0
%         n=n-1;
%     end
%     x1=interp1([y(n) y(n+1)],[x(n) x(n+1)],maxy/2);
%     slope1=(y(n+1)-y(n))./(x(n+1)-x(n));
%     
%     n=indmax;
%     while oy(n)>0
%         n=n+1;
%     end
%     x2= interp1([y(n-1) y(n)],[x(n-1) x(n)],maxy/2);
%     slope2=(y(n-1)-y(n))./(x(n-1)-x(n));
%     
%     FWHM=x2-x1;
%     hwhm1=xo-x1;
%     hwhm2=x2-xo;
% catch
%     FWHM=NaN;
% end
% 
% function [index,closestval]=val2ind(x,val)
% % Returns the index and the value of the element of vector x that is closest to val
% % If more than one element is equally close, returns vectors of indicies and values
% % Tom O'Haver (toh@umd.edu) October 2006
% % Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% % [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
% dif=abs(x-val);
% index=find((dif-min(dif))==0);
% closestval=x(index);
% function Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth)
% % Resolution enhancement function by derivative method. the
% % arguments factor1 and factor 2 are 2nd and 4th derivative weighting
% % factors. Larger values of factor1 and factor2 will reduce the 
% % peak width but will cause artifacts in the baseline near 
% % the peak.  Adjust the factors for the the best compromise. 
% % Use minimum smooth width needed to reduce excess noise. 
% % Use iSignal to adjust these factors 
% % 
% % Example:
% %  x=[0:.01:18];
% %  y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2);
% %  Enhancedsignal=enhance(y,1000,1000000,3)
% %  plot(x,y,x,Enhancedsignal,'r')
% %
% % Related functions: SegmentedSharpen.m
% 
% d2=secderiv(signal);  % Computes second derivative
% d4=secderiv(d2);   % Computes fourth derivative
% Enhancedsignal = signal-factor1.*tsmooth(d2,SmoothWidth)+...
% factor2.*tsmooth(tsmooth(tsmooth(d4,SmoothWidth),SmoothWidth),SmoothWidth);
% 
% function d=secderiv(a)
% % Second derivative of vector using 3-point central difference.
% % Example: secderiv([0 1 2 3 4 3 2 1 0]) yields [ 0 0 0 0 -2 0 0 0 0]
% %  T. C. O'Haver, 2006.
% n=length(a);
% d=zeros(size(a));
% for j = 2:n-1
%   d(j)=a(j+1) - 2.*a(j) + a(j-1);
% end
% d(1)=d(2);
% d(n)=d(n-1);
% 
% function s=tsmooth(Y,w)
% %  tsmooth(Y,w) smooths vector Y by a triangular function of halfwidth w
% %  T. C. O'Haver, 1988.
% v=ones(1,w);v=conv(v,v);
% S=conv(Y,v);
% startpoint=(length(v) + 1)/2;
% endpoint=length(Y)+startpoint-1;
% 
% s=S(startpoint:endpoint) ./ sum(v);
% 
% % Copyright (c) 2018, Thomas C. O'Haver
% % 
% % Permission is hereby granted, free of charge, to any person obtaining a copy
% % of this software and associated documentation files (the "Software"), to deal
% % in the Software without restriction, including without limitation the rights
% % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % copies of the Software, and to permit persons to whom the Software is
% % furnished to do so, subject to the following conditions:
% % 
% % The above copyright notice and this permission notice shall be included in
% % all copies or substantial portions of the Software.
% % 
% % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% % THE SOFTWARE.