function DerivativeNumericalPrecisionDemo
% Self-contained function that shows how the numerical precision limits of
% the computer effects the first through fourth derivatives of a very
% finely sampled ("noiseless") Gaussian band, showing both the waveforms
% (figure(1) and their freqency spectra (Figure(2). The numerical precision
% limits of the computer create random noise at very high frequencies,
% which are emphasized by differentiation, and by fourth derivative that
% noise overwhelms the signal frequencies at lower frequencies. Smoothing
% with three passes of a sliding average smooth with a smooth ratio of 0.2
% removes most of the noise.

% Create an independent variable ("x-axis")
increment=.003; % A small x-axis increment pushes the numerical precision limits
x=1:increment:300;

% Create a Gaussian ("bell-shaped") signal as a function of x, centered at
% x=150 with width 50.
y=gaussian(x,150,50);

% Plot the 1st, 2nd, 3rd, and 4th derivatives in the four quadrants of the
% plot window
figure(1)
subplot(2,2,1)
dy=deriv(y)./increment;
plot(x,dy)
title('First Derivative')
subplot(2,2,2)
dy2=deriv2(y)./(increment.^2);
plot(x,dy2)
title('Second Derivative')
subplot(2,2,3)
dy3=deriv3(y)./(increment.^3);
plot(x,dy3)
title('Third Derivative')
subplot(2,2,4)
dy4=deriv4(y)./(increment.^4);
plot(x,dy4)
title('Fourth Derivative')

% Prepare log-log plot of the frequency spectra of the four derivatives in
% the four quadrants of the plot window
figure(2);
clf
subplot(2,2,1)
PlotFrequencySpectrum(x,dy',4,0,0);
title('Frequency spectrum of First Derivative')
subplot(2,2,2)
PlotFrequencySpectrum(x,dy2',4,0,0);
title('Frequency spectrum of Second Derivative')
subplot(2,2,3)
PlotFrequencySpectrum(x,dy3',4,0,0);
title('Frequency spectrum of Third Derivative')
subplot(2,2,4)
PlotFrequencySpectrum(x,dy4',4,0,0);
title('Frequency spectrum of Fourth Derivative')

figure(3)
hold on;PlotFrequencySpectrum(x,dy',4,0,0);
PlotFrequencySpectrum(x,dy2',4,0,0);
PlotFrequencySpectrum(x,dy3',4,0,0);
PlotFrequencySpectrum(x,dy4',4,0,0);
hold off

% -----------------------------------------------------
% Internal sub-functions

function d=deriv(a)
% First derivative of vector using 2-point central difference.
% Example: deriv([1 1 1 2 3 4]) yeilds [0 0 .5 1 1 1]
%  T. C. O'Haver, 1988.
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j));
end

function d2=deriv2(a)
% Second derivative of vector using 3-point central difference.
%  T. C. O'Haver, 2006.
n=length(a);
d2=zeros(size(a));
for j = 2:n-1;
  d2(j)=a(j+1) - 2.*a(j) + a(j-1);
end
d2(1)=d2(2);
d2(n)=d2(n-1);

function d3=deriv3(a)
% Third derivative of vector a
%  T. C. O'Haver, 2008.
n=length(a);
d3=zeros(size(a));
for j = 3:n-2;
  d3(j)=a(j+2) - 2.*a(j+1) + 2.*a(j-1) - a(j-2);
end
d3(1:2)=d3(3);
d3(n-1:n)=d3(n-2);

function d4=deriv4(a)
% 4th derivative of vector a
%  T. C. O'Haver, 2008.
% Example: Shows numerical round-off error in Matlab; a smaller increment
% than 0.0003 causes increasing noise. 
% x=-5:.0003:5; y=exp(-(x).^2);plot(deriv4(y))
n=length(a);
d4=zeros(size(a));
for j = 3:n-2;
  d4(j)=a(j+2) - 4.*a(j+1) + 6.*a(j) - 4.*a(j-1) + a(j-2);
end
d4(1:2)=d4(3);
d4(n-1:n)=d4(n-2);

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);

function PowerSpectrum=PlotFrequencySpectrum(xx,yy,plotmode,XMODE,LabelPeaks)
% Computes the power spectrum of x,y with optional linear/log plotting.
% plotmode: =0: no plot; =1:linear, =2:semilog X, =3:semilog Y; =4: log-log)
% XMODE: =0 for frequency Spectrum (x is frequency); =1 for periodogram (x is time).
% If LabelPeaks = 1, uses findpeaksx to locate and label frequency of peaks
% Change peak detection parameters in lines 51-54.
% (c) Tom O'Haver 2014.
% EXAMPLE:
% x=[0:.01:2*pi]';
% y=sin(200*x)+randn(size(x));
% subplot(2,1,1);
% plot(x,y);
% subplot(2,1,2);
% PowerSpectrum=PlotFrequencySpectrum(x,y,1,0,1);
fy=fft(yy);
sy=fy .* conj(fy); % Compute power spectrum
plotrange=1:length(fy)/2;
if XMODE,
    f=range(xx)./(plotrange-1);
else
    f=((plotrange-1)./range(xx));
end
realsy=real(sy(plotrange));
maxpower=max(realsy);
maxx=val2ind(realsy,maxpower);
PowerSpectrum=([f' realsy]);
if plotmode,
    
    switch plotmode,
        case 1
            plot(f,realsy,'r.-')
            ylabel('Linear y')
        case 2
            semilogx(f,realsy,'r.-')
            ylabel('Linear y')
        case 3
            semilogy(f,realsy,'r.-')
            ylabel('Log y')
        case 4
            loglog(f,realsy,'r.-')
            ylabel('Log y')
        otherwise,
    end
    if XMODE,
        title('Periodogram: x=Time. ')
        xlabel('Time')
    else
        title('Spectrum: x=Frequency (e.g. 1/time)')
        xlabel('Frequency')
    end
    if LabelPeaks,
        AmpT=0.05.*max(PowerSpectrum(:,2));
        SlopeT=0.01;
        SmoothW=3;
        FitW=3;
        P=findpeaksx(PowerSpectrum(:,1),PowerSpectrum(:,2),SlopeT,AmpT,SmoothW,FitW);
        text(P(:,2),P(:,3),num2str(P(:,2)))
    end
end

function P=findpeaksx(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% function P=findpeaksx(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Function to locate the positive peaks in a noisy x-y time series data
% set. Detects peaks by looking for downward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (P) containing peak number and position and height of each peak. Arguments "slopeThreshold",
% "ampThreshold" and "smoothwidth" control peak sensitivity.
% Higher values will neglect smaller features. "Smoothwidth" is
% the width of the smooth applied before peak detection; larger
% values ignore narrow peaks. If smoothwidth=0, no smoothing
% is performed. "Peakgroup" is the number points around the top
% part of the peak that are taken for measurement. The 
% the local maximum is taken as the peak height and position.
% The argument "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar)
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver.  Version 1, April, 2016  
%
% Examples:
% findpeaksx(0:.01:2,humps(0:.01:2),0,-1,5,5)
% x=[0:.01:50];findpeaksx(x,cos(x),0,-1,5,5)
% x=[0:.01:5]';y=x.*sin(x.^2).^2;plot(x,y);findpeaksx(x,y,0,-1,5,5)
% x=[-10:.1:10];y=exp(-(x).^2);findpeaksx(x,y,0.005,0.3,3,5,3);
%
% Find 12000 peaks in less than one second
% x=[0:.01:500]';
% y=x.*sin(x.^2).^2;
% tic;
% P=findpeaksx(x,y,0,0,3,3);
% toc;
% NumPeaks=max(P(:,1))
%
% Find, plot, and number noisy peaks with unknown positions
% x=-50:.2:50;
% y=exp(-(x).^2)+exp(-(x+50*rand()).^2)+.02.*randn(size(x));
% plot(x,y,'m.')
% P=findpeaksx(x,y,0.001,0.2,1,2,3);
% text(P(:,2),P(:,3),num2str(P(:,1)))
% disp('           peak #    Position    Height')
% disp(P)
%
% Related functions:
% findvalleys.m, findpeaksG.m, findpeaksb.m, findpeaksb3.m,
% findpeaksplot.m, peakstats.m, findpeaksnr.m, findpeaksGSS.m,
% findpeaksLSS.m, findpeaksfit.m, findsteps.m, findsquarepulse.m, idpeaks.m

% Copyright (c) 2013, 2014 Thomas C. O'Haver
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
%
if nargin~=7;smoothtype=1;end  % smoothtype=1 if not specified in argument
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end
if smoothwidth<1;smoothwidth=1;end
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=fastsmooth(deriv(y),smoothwidth,smoothtype);
else
    d=deriv(y);
end
n=round(peakgroup/2+1);
P=[0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
PeakX=0;
PeakY=0;
% plot(x,d,'.')
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth-1,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold, % if slope of derivative is larger than SlopeThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup<3;
                    PeakY=max(yy);
                else
                    PeakY=mean(yy);
                end
                pindex=val2ind(yy,PeakY);
                PeakX=xx(pindex(1));
            % Construct matrix P. One row for each peak detected,
            % containing the peak number, peak position (x-value) and
            % peak height (y-value). If peak measurement fails and
            % results in NaN, or if the measured peak height is less
            % than AmpThreshold, skip this peak
            if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold,
                % Skip this peak
            else % Otherwise count this as a valid peak
                P(peak,:) = [round(peak) PeakX PeakY];
                peak=peak+1; % Move on to next peak
            end
        end
    end
end
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth
%  of width w. Version 2.0, May 2008.
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar)
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively
%     smaller smooths the closer to the end. (In this mode the
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%  T. C. O'Haver, May, 2008.
if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
switch type
    case 1
        SmoothY=sa(Y,w,ends);
    case 2
        SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
        SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
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
% ----------------------------------------------------------------------
