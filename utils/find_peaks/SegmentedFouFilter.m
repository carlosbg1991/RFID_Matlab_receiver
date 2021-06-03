function ffSignal = SegmentedFouFilter(y,samplingtime,centerFrequency,filterWidth,shape,mode )
% ffSignal = SegmentedFouFilter(y,samplingtime,centerFrequency,filterWidth,shape,mode)
% Segmented Fourier filter function for time-series signal vector y; 'samplingtime'
% is the total duration of sampled signal in sec, millisec, or microsec;
% 'centerfrequency' and 'filterWidth' are the center frequency and width 
% of the filter in Hz, KHz, or MHz, respectively, and must be vectors of equal length. 'Shape' determines the 
% sharpness of the cut-off. If shape = 1, the filter is Gaussian; as  
% shape increases the filter shape becomes more and more rectangular. 
% Set mode = 0 for band-pass filter, mode = 1 for band-reject (notch) filter.  
% FouFilter returns the filtered signal.
%
% Example 1: Variable frequency bandpass filter with f=10 sinewave input.
% x=0:.001:2*pi;
% f=10;
% y=sin(2*pi*f*x);
% sff=SegmentedFouFilter(y,2*pi,[8 9 10 11 12],[1 1 1 1 1],3,0);
% plot(x,sff)
%
% Example 2: low-pass filter with decreasing bandwidth
% load WideningPeaks
% sff=SegmentedFouFilter(y,max(x),[0 0 0 0 0],[.03 .025 .02 .015 .01],3,0);
% plot(x,y,'c',x,sff,'k')
% legend('Original signal','Fourier filtered')
%
% Example 3: As example 2, with gradient bandwidth
% load WideningPeaks
% NumSegments=10;
% endv=10;startv=30;
% vstep=(endv-startv)/NumSegments;
% centerFrequencies=startv:vstep:endv;
% shape=3;
% sff=SegmentedFouFilter(y,10,ones(1,NumSegments),centerFrequencies,shape,0);
% plot(x,z,x,sff)
% legend('Noiseless data','Filtered noisy data')
%
%  T. C. O'Haver (toh@umd.edu), November 2019. Based on version 2 of
%  FouFilter.m, which included a correction thanks to Dr. Raphael Attié,
%  NASA/Goddard Space Flight Center
%   Copyright (c) 2019, Thomas C. O'Haver
   
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
NumSegments=length(centerFrequency);
SegLength=round(ly./NumSegments);
ffSegment=zeros(ly,NumSegments);
ffSignal=zeros(1,ly);
for Segment=1:NumSegments
    ffSegment(:,Segment)=ff(y,samplingtime,centerFrequency(Segment),filterWidth(Segment),shape,mode);
    startindex=(1+(Segment-1)*SegLength);
    endindix=startindex+SegLength-1;
    if endindix>ly,endindix=ly;end
    indexrange=startindex:endindix;
    ffSignal(indexrange)=ffSegment(indexrange,Segment);
end

function ry=ff(y,samplingtime,centerFrequency,filterWidth,shape,mode)
center=centerFrequency*samplingtime; %  center harmonic (fourier component)
width=filterWidth*samplingtime; %  width of filter (in harmonics)

fy=fft(y); % Fourier transform of signal

% lft1=1:(length(fy)/2);
% lft2=(length(fy)/2+1):length(fy);
% 
% ffilter1=ngaussian(lft1,center+1,width,shape);
% ffilter2=ngaussian(lft2,length(fy)-center+1,width,shape);

lft1=0:floor(length(fy)/2)-1;
lft2=floor(length(fy)/2):length(fy)-1;

ffilter1=ngaussian(lft1,center,width,shape);
ffilter2=ngaussian(lft2,length(fy)-center,width,shape);

ffilter=[ffilter1,ffilter2];
if mode==1
    ffilter=1-ffilter; 
end
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter;  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Recover filter signal from Fourier transform

function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = flattened Gaussian centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
% Shape is Gaussian when n=1. Becomes more rectangular as n increases.
%  T. C. O'Haver, 1988, revised 2014
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
g = exp(-((x-pos)./(0.6005612.*wid)) .^(2*n));