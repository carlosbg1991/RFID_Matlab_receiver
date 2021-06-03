function [xx,yy,irange]=plotxrange(x,y,x1,x2)
% function [xx,yy,irange]=plotxrange(x,y,x1,x2)
% Extracts and plots values of vectors x,y only for x values between x1 and
% x2. Returns extracted values in vectors xx,yy and the range of index
% values in irange. Ignores values of x1 and x2 outside the range of x.
%
% Example:
% x=1:.1:100;
% y=x.*sin(x).^2;
% plotxrange(x,y,5,20);
%
% Related functions: segplot.m, val2ind.m; plotit.m

% Copyright (c) 2019, Thomas C. O'Haver
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

irange=val2ind(x,x1):val2ind(x,x2);
xx=x(irange);
yy=y(irange); 
plot(xx,yy)

function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is
% closest to val. If more than one element is equally close, returns vectors
% of indicies and values.
% Tom O'Haver (toh@umd.edu) October 2006
%
% Example 1: 
% If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
%
% Example 2:
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
%
% Example 3:
% Find the value of x between 0 and 6 for which sin(x) is maximum.
% x=0:.001:6;
% y=sin(x);
% x0=x(val2ind(y,max(y)))
% (returns 1.5708, which as expected is pi/2).

dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);