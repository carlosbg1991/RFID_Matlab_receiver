function powertest2
% Test of power method for measuring area of small shouldering peak
% You can change parameters in lines 6-15
format compact; format short g
clear
increment=.05;
x=1:increment:50;
HL=4; % height of large peak to the left
Ho=.3; % peak height of measured (second) peak
p1=34; % position of first peak
p2=3; % position of measured (second) peak
w1=3.1; % width of first peak
w2=3.1; % width of measured (second) peak
baseline=.01; % Baseline as fraction of height of second peak
power=22;
%
x1=val2ind(x,p1); % Index value of peak position 1
x2=val2ind(x,p2); % Index value of peak position 2
g1=HL.*gaussian(x,p1,w1); % gaussian, height HL, centered at p1, FWHM=w1
g2=Ho.*gaussian(x,p2,w2); % gaussian, height Ho, centered at p2, FWHM=w2
y=g1+g2+baseline.*Ho; % Sum of the two peaks plus baseline
ValleyPoint=min(y(x1:x2));  % Y value of valley between peaks
ValleyIndex=val2ind(y,ValleyPoint); % Index of valley
StartPoint=ValleyIndex; % Start point to computer perp drop areas
EndPoint=val2ind(x,45); % Last point in signal
TrueArea=sum(g2)*increment; % Actual area under peak 2 without peak 1
PerpDropArea=sum(y(StartPoint:EndPoint))*increment; % Perp drop area of peak 2
% Power computations
yp=y.^power; % yp = y raised to power
ValleyMinY=min(yp(x1:x2)); % Y^2 value of valley between peaks
ValleyIndex=val2ind(yp,ValleyMinY); % Index of valley y^2
An=sum(yp(StartPoint:EndPoint))*increment; % Perp drop area of peak 2 ^ power
nyp=yp/max(yp(StartPoint:EndPoint)); % Normzalized y^power
nAn=sum(nyp(StartPoint:EndPoint)).*increment; % Norm Perp drop area
PerpDropAreaPowerNorm = Ho*nAn*sqrt(power); % Equation 5 in paper
PercentErrorNorm=100.*(PerpDropArea-TrueArea)./TrueArea;
PercentErrorPower=100.*(PerpDropAreaPowerNorm-TrueArea)./TrueArea;
clf;
plot(x,y./max(y),'g--',x,yp./max(yp(x2:length(x))),'r--',x(StartPoint:EndPoint),y(StartPoint:EndPoint)./max(y),'-g')
axis([25 50 0 1])
hyp=halfwidth(x,yp); % halfwidth of yp
ValleyDip=(Ho-ValleyPoint)/Ho; % fractional dip in valley between peaks
hold on;
% plot([x(ValleyIndex) x(ValleyIndex)],[0 1]);
plot(x(val2ind(x,p2-2*hyp):val2ind(x,p2+2*hyp)),nyp(val2ind(x,p2-2*hyp):val2ind(x,p2+2*hyp)),'-r');
hold off
xlabel('X or time.   Solid lines are regions of area measurement for second peak.')
ylabel('Normalized y or y^p')
title('Green : y (Normalized to peak 1)      Red : y^p (Normalized to peak 2)')
disp(' ')
disp('    Height ratio  Valley dip     Baseline       Power'); 
disp([HL/Ho ValleyDip baseline power])
disp('% error of small peak area measurement for second peak.')
disp('       Normal       Power')
disp([PercentErrorNorm PercentErrorPower])

%-------- functions----------------------------------
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

