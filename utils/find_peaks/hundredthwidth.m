function [FWHM,slope1,slope2]=hundredthwidth(x,y,xo)
% function FWHM=hundredthwidth(x,y,xo) computes the full width at 1/10 maximum
% of the maximum of any shape peak that peaks near xo and has a zero
% baseline. If xo is omitted, it computes the halfwidth from the maximum y
% value. Not highly accurate if the function is too sparsely sampled. Also
% returns the leading and trailing edge slopes. Tom O'Haver (toh@umd.edu)
% Version 1, August 2016
%
% Example 1:
% x=-5:.01:5;
% y=sinc(x).^2;
% FWHM=hundredthwidth(x,y,0) % Area of central peak
% FWHM=hundredthwidth(x,y,1.5) % Area of positive side peak
%
% Example 2:
% x=[0:.1:10];
% W=3; % W is the true half-width
% y=gaussian(x,5,W);
% FWHM=hundredthwidth(x,y,5)
% FWHM=hundredthwidth(x,y) % gives the same result since there is only one peak
%
if nargin==2,xo=x(val2ind(y,max(y)));end
try   
    indmax=val2ind(x,xo);
    maxy=y(indmax);
    oy=y-maxy/100;
    
    n=indmax;
    while oy(n)>0,
        n=n-1;
    end
    x1=interp1([y(n) y(n+1)],[x(n) x(n+1)],maxy/100);
    slope1=(y(n+1)-y(n))./(x(n+1)-x(n));
    n=indmax;
    while oy(n)>0,
        n=n+1;
    end
    x2= interp1([y(n-1) y(n)],[x(n-1) x(n)],maxy/100);
    slope2=(y(n-1)-y(n))./(x(n-1)-x(n));
    FWHM=x2-x1;
    
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

% Copyright (c) 2012, Thomas C. O'Haver
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