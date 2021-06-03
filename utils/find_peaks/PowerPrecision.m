function PowerPrecision
% Test of power method for measuring area of small shouldering peak
% You can change parameters in lines 6-15
format compact; format short g
clear
increment=.05; % difference bwtween adjacent x-axis values
x=1:increment:50; % independent variable axis (e.g. time)
Ho=0.2; % peak height of measured (second) peak
HL=4; % height of large peak to the left
p1=34; % position of first peak
p2=40; % position of measured (second) peak
w1=3; % width (FWHM) of larger (first) peak
w2=3; % width (FWHM) of measured (second) peak
g1=HL.*gaussian(x,p1,w1); % gaussian, height HL, centered at p1, FWHM=w1
g2=Ho.*gaussian(x,p2,w2); % gaussian, height Ho, centered at p2, FWHM=w2
baseline=.05; % Baseline as fraction of height of second peak
power=8; % Power to which the data are raised
StartPlot=25; % X value to start plotting
noise=.001; % Standard deviation of random white noise added to signal
SmoothWidth=11; % Number of data points in smooth
w1s= halfwidth(x,fastsmooth(g1,SmoothWidth)); % width (FWHM) after smoothing
w2s= halfwidth(x,fastsmooth(g2,SmoothWidth)); % width (FWHM) after smoothing
Rsmoothed=1.18.*(p2-p1)./(w1s+w2s); % Resolution after smoothing
disp(' ') % Print out settings in command window
disp('   Height ratio    Resolution     Baseline')
disp([HL/Ho Rsmoothed baseline])
disp('   Noise   SmoothWidth      Power')
disp([noise SmoothWidth power])
for n=1:20
    x1=val2ind(x,p1); % Index value of peak position 1
    x2=val2ind(x,p2); % Index value of peak position 2
    y=g1+g2+baseline.*Ho+noise.*randn(size(x)); % Sum of the two peaks plus baseline
    if SmoothWidth>0,y=fastsmooth(y,SmoothWidth);end
    ValleyPoint=min(y(x1:x2));  % Y value of valley between peaks
    ValleyIndex=val2ind(y,ValleyPoint); % Index of valley
    StartPoint=ValleyIndex; % Start point to computer perp drop areas
    EndPoint=val2ind(x,45); % Last point in signal
    TrueArea=sum(g2)*increment; % Actual area under peak 2 without peak 1
    PerpDropArea(n)=sum(y(StartPoint:EndPoint))*increment; % Perp drop area of peak 2
    % Power computations
    yp=y.^power; % yp = y raised to power
    ValleyMinY=min(yp(x1:x2)); % Y^2 value of valley between peaks
    ValleyIndex=val2ind(yp,ValleyMinY); % Index of valley y^2
    An=sum(yp(StartPoint:EndPoint))*increment; % Perp drop area of peak 2 ^ power
    nyp=yp/max(yp(StartPoint:EndPoint)); % Normzalized y^power
    nAn=sum(nyp(StartPoint:EndPoint)).*increment; % Norm Perp drop area
    PerpDropAreaPowerNorm(n) = Ho*nAn*sqrt(power); % Equation 5 in paper
    % Plot the peaks in the figure window
    clf;
    plot(x,y./max(y),'g--',x,yp./max(yp(x2:length(x))),'r--',x(StartPoint:EndPoint),y(StartPoint:EndPoint)./max(y),'-g')
    axis([StartPlot 50 0 1]);
    hyp=halfwidth(x,yp); % halfwidth of yp
    ValleyDip=(Ho-ValleyPoint)/Ho; % fractional dip in valley between peaks
    hold on;
    % plot([x(ValleyIndex) x(ValleyIndex)],[0 1]);
    plot(x(val2ind(x,p2-2*hyp):val2ind(x,p2+2*hyp)),nyp(val2ind(x,p2-2*hyp):val2ind(x,p2+2*hyp)),'-r');
    hold off
    xlabel('X or time.   Solid lines are regions of area measurement for second peak')
    ylabel('Normalized y or y^p')
    title('Green : y (Normalized to peak 1)      Red : y^p (Normalized to peak 2)')
    drawnow;
end
MeanStandard=mean(PerpDropArea);
RSDstandard=100.*(std(PerpDropArea)./MeanStandard);
MeanPower=mean(PerpDropAreaPowerNorm);
RSDpower=100.*(std(PerpDropAreaPowerNorm)./MeanPower);
PercentErrorNorm=100.*(MeanStandard-TrueArea)./TrueArea;
PercentErrorPower=100.*(MeanPower-TrueArea)./TrueArea;
% Print out results in command window
disp('Percent error')
disp('       Normal       Power')
disp([PercentErrorNorm PercentErrorPower])
disp('Percent Relative standard deviation')
    disp('       Normal       Power')
    disp([RSDstandard RSDpower])

% Peakfit.m applied to raw y data, for comparison
% [FitResults,GOF]=peakfit([x;y],(p1+p2)/2,2*(p2-p1+w1+w2),2,1,0,10,0,1,0,0)
% PFArea=FitResults(2,5) % Peakfit area of peak 2
% PercentErrorPeakfit=100.*(PFArea-TrueArea)./TrueArea



%-------- subfunctions----------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);

function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is
% closest to val If more than one element is equally close, returns vectors
% of indicies and values.
% Tom O'Haver (toh@umd.edu) October 2006
%
% Example 1:
% If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
%
% Example 2:
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]

% Copyright (c) 2012, Thomas C. O'Haver
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

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
    while oy(n)>0,
        n=n-1;
    end
    x1=interp1([y(n) y(n+1)],[x(n) x(n+1)],maxy/2);
    slope1=(y(n+1)-y(n))./(x(n+1)-x(n));
    
    n=indmax;
    while oy(n)>0,
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
