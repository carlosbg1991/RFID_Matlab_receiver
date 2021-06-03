function [s,xx,yy]=segplot(x,y,NumSegs,seg)
% [s,xx,yy]=segplot(x,y,NumSegs,seg) divides y into "NumSegs" equal-length
% segments and plots the x, y data with segments marked by vertical lines,
% each labeled with a small segment number at the bottom. Returns a vector
% 's' of segment indexes, and the subset xx,yy, of values in the segment
% number 'seg'. If the 4th input argument, 'seg', is included, it plots
% this particular segment only.
%
% Example:  
% x=1:.1:100;
% y=x.*sin(x).^4;
% figure(1)
% [s,xx,yy]=segplot(x,y,41);
% figure(2)
% [s,xx,yy]=segplot(x,y,41,10);
%
% Related functions: plotxrange.m, val2ind.m, SegmentedSmooth.m
 
%   Copyright (c) 2017, Thomas C. O'Haver
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, includingnwithout limitation the rights
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
SegLength=round(ly./NumSegs);
s=zeros(1,NumSegs);
if nargin==3
    plot(x,y,'k')
    hold on
    for Segment=1:NumSegs
        startindex=(1+(Segment-1)*SegLength);
        endindix=startindex+SegLength-1;
        if endindix>ly,endindix=ly;end
        plot([x(endindix) x(endindix)],[min(y) max(y)],'c')
        s(Segment)=startindex;
        text(x(startindex),min(y),num2str(Segment),'FontSize',6 )
    end
    hold off
    xx=0;
    yy=0;
end
if nargin==4
    for Segment=1:NumSegs
        startindex=(1+(Segment-1)*SegLength);
        s(Segment)=startindex;
    end
    xx=x(s(seg):s(seg+1));
    yy=y(s(seg):s(seg+1));
    plot(xx,yy)
    title(['Segment number ' num2str(seg)])
end


