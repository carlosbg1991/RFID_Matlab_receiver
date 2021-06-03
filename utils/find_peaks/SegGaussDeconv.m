function Deconvoluted=SegGaussDeconv(x,y,w)
% SegGaussDeconvPlot(x,y,w) divides y into a number of equal-length
% segments defined by the length of the vector 'w', then each segment is
% deconvoluted with a Gaussian where t is corresponing element of the
% vector w. Any number and sequence of t values can be used.
% Related functions:
% SegExpDeconv, SegExpDeconvPlot, SegGaussDeconv, SegGaussDeconvPlot
 
%   Copyright (c) 2017, Thomas C. O'Haver
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
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
NumSegments=length(w);
SegLength=round(ly./NumSegments);
DeconvSegment=zeros(ly,NumSegments);
Deconvoluted=zeros(1,ly);
% figure(8);clg
for Segment=1:NumSegments,
    lx=length(x);hx=round(lx/2);
    c=[gaussian(x(1:hx),0,w(Segment));gaussian(x(hx+1:lx),max(x),w(Segment))];   % Gaussian convolution function, c
    DeconvSegment(:,Segment)=ifft(fft(y)./fft(c)).*sum(c); % DeconvSegment = y deconvoluted
    startindex=(1+(Segment-1)*SegLength);
    endindix=startindex+SegLength-1;
    if endindix>ly,endindix=ly;end
    indexrange=startindex:endindix;
    Deconvoluted(indexrange)=DeconvSegment(indexrange,Segment);
end


