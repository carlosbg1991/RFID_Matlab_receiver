function DerivativeDemo 
% DerivativeDemo.m is a self-contained Matlab/Octave function that uses the
% ProcessSignal and plotit.m to demonstrate an application of
% differentiation to quantitative analysis of a peak buried in an unstable
% background.  To run it, just type DerivativeDemo at the command prompt.
% You can change several of the internal variables (e.g. Noise, in
% line 12 or BackgroundAmplitude in line 13), etc. 
clf
xincrement=.1;
x=1:xincrement:256;
SignalWidth=30;
Noise=.002;
BackgroundAmplitude=.03;
BackgroundPeak=-100;
BackgroundWidth=600;
SignalAmplitudes=[0 .005 .01 .02 .03 .04];
DerivativeMode=2;
SmoothMode=3;
SmoothWidth=201;
SmoothRatio=201/(SignalWidth/xincrement);
em=0;

y1=rand(1)*BackgroundAmplitude.*gaussian(x,BackgroundPeak,BackgroundWidth)+SignalAmplitudes(1).*gaussian(x,128,SignalWidth)+Noise.*randn(size(x));
y2=rand(1)*BackgroundAmplitude.*gaussian(x,BackgroundPeak,BackgroundWidth)+SignalAmplitudes(2).*gaussian(x,128,SignalWidth)+Noise.*randn(size(x));
y3=rand(1)*BackgroundAmplitude.*gaussian(x,BackgroundPeak,BackgroundWidth)+SignalAmplitudes(3).*gaussian(x,128,SignalWidth)+Noise.*randn(size(x));
y4=rand(1)*BackgroundAmplitude.*gaussian(x,BackgroundPeak,BackgroundWidth)+SignalAmplitudes(4).*gaussian(x,128,SignalWidth)+Noise.*randn(size(x));
y5=rand(1)*BackgroundAmplitude.*gaussian(x,BackgroundPeak,BackgroundWidth)+SignalAmplitudes(5).*gaussian(x,128,SignalWidth)+Noise.*randn(size(x));
y6=rand(1)*BackgroundAmplitude.*gaussian(x,BackgroundPeak,BackgroundWidth)+SignalAmplitudes(6).*gaussian(x,128,SignalWidth)+Noise.*randn(size(x));
y=[y1;y2;y3;y4;y5;y6];
plotrange=2*SmoothWidth:length(x)-2*SmoothWidth;
PeakPoint=128./xincrement;
figure(1)
subplot(2,2,1);plot(x(plotrange),y(:,plotrange))
title('Original Signals');
for n=1:6,
   sy(n,:)=ProcessSignal(x,y(n,:),0,SmoothWidth,SmoothMode,em,0,0,0,0,0);
   dy(n,:)=ProcessSignal(x,y(n,:),DerivativeMode,SmoothWidth,SmoothMode,em,0,0,0,0,0);
   da(n)=-dy(n,PeakPoint);
end
subplot(2,2,2);plot(x(plotrange),dy(:,plotrange))
title('Smoothed Derivative Signals');

subplot(2,2,3)
plotit(SignalAmplitudes,sy(:,PeakPoint)',1);
title('Amplitude of smoothed signal at center of peak');
xlabel('Actual peak amplitude')
ylabel('Measured peak amplitude')

subplot(2,2,4)
plotit(SignalAmplitudes,da,1);
title('Derivative results');
xlabel('Actual peak amplitude')
ylabel('Measured derivative peak amplitude')
%-----------------------------------------------------------------
function [coef,RSquared,StdDevs,BootResults]=plotit(xi,yi,polyorder,datastyle,fitstyle)
% [coef,RSquared,StdDevs]=plotit(x,y,n,datastyle,fitstyle)
% plotit accepts data in the form of a single vector "y", a pair of vectors
% "x" and "y", or a 2xn or nx2 matrix with x in first row or column. If the
% input argument "polyorder" is supplied, it fits the data to a polynomial
% of order "polyorder", plots the data in red dots and the fit as a blue
% line, and displays the fit coefficients and R-squared in the upper left
% corner of the graph. Polyorder=1 for straight line, =2 for quadratic
% (parabola) etc. 
%   The 4th and 5th input arguments are strings defining the line and symbol
% style and color, in standard Matlab convention. The strings, in single
% quotes, are made from one element from any or all the following 3
% columns:
%  
%            b     blue          .     point              -     solid
%            g     green         o     circle             :     dotted
%            r     red           x     x-mark             -.    dashdot 
%            c     cyan          +     plus               --    dashed   
%            m     magenta       *     star             (none)  no line
%            y     yellow        s     square
%            k     black         d     diamond
%            w     white         v     triangle (down)
%                                ^     triangle (up)
%                                <     triangle (left)
%                                >     triangle (right)
%                                p     pentagram
%                                h     hexagram
% For example, plotit(x,y,3,'or','-g') plots the data as red circles and
% the fit as a green solid line.
%   If the 4th output argument (BootResults) is supplied, plotit computes
% coefficient error estimates by the bootstrap method and returns the
% results in the matrix "BootResults" (of size 5 x polyorder+1). You can
% change the number of bootstrap samples in line 93. 
% Tom O'Haver, Version 6, January 2015, added optional specification of
% plot style and color. Version 5 added standard deviations of coefficients
% as 3rd output argument (moved bootstrap results to 4th output argument).
% Examples:
% x=[1 2 3 4 5];y=[0 2 3 3 5];
% plotit(y); % plot y only vs index number, no fit.
% plotit(x,y); % plot x vs y only, data in separate vectors
% plotit([x;y]); % plot data only, data in matrix with x in first row
% plotit([x;y]'); % plot data only, data in matrix with x in first column
% plotit(y,2); % plot y vs index and fit to second order polynomial
% plotit(x,y,2); % plot x vs y and fit t  second order polynomial
% plotit(x,y,3,'or'); % as above, plotting data as red circles
% plotit(x,y,3,'or','-g'); % data as red circles, fit as green line
% plotit([x;y],3);  % plot and fit data in matrix
% [coef, RSquared]=plotit([x;y],2) % return fit coefficients and r-squared
% [coef,RSquared,STD]=plotit([x;y],2) Returns vector standard devations
% of coefficients as 'STD'. (STD./coef computes relative standard deviation)
% [coef, RSquared,,StdDevs,BootResults]=plotit(x,y,polyorder) % computes
% coefficient error estimates by the bootstrap method and returns the
% results in the matrix "BootResults" 
%
% You can use plotit.m to linearize and plot nonlinear relationships, such as
% y = a*exp(b*x) ; [coeff,R2]=plotit(x,log(y),1);a=exp(coeff(2));b=coeff(1);
% y = a*ln(b*x) ; [coeff,R2]=plotit(log(x),y,1);a=coeff(1);b=log(coeff(2));
% y=ax^b ; [coeff,R2]=plotit(log(x),log(y),1);a=exp(coeff(2));b=coeff(1);
%
% Related functions (http://tinyurl.com/cey8rwh): 
% plotfita.m animates the bootstrap process for instructional purposes. 
% logplotfit.m plots and fits log(x) vs log(y), for data that follows a
%    power law relationship or that covers a very wide numerical range. 
% trypoly(x,y) fits a series of polynomials to the data x,y, for
%   polynomial orders 1 through length(x)-1, returns the R2 values in a
%   vector. Shows that, for any data, R2 -> 1 as order -> length(x)-1.
% trydatatrans(x,y,order) tries 8 different simple data transformations
%   on the data x,y, fits the transformed data to a polynomial of order
%   'polyorder', returns the R2 values in a vector. Displays results in
%   3 x 3 array of small plots.

% Copyright (c) 2015, Thomas C. O'Haver

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
% process input arguments
NumTrials=100;
if nargin==1, % Single argument, might be y vector or x,y matrix
    polyorder=0;
    datasize=size(xi);
    if isvector(xi), % xi is vector, not matrix
        y=xi;
        x=1:length(y); % Use this only to create an x vector if needed
        xlabelstring='index number';
    else % x is a matrix, split it into x and y vectors
        if datasize(1)<datasize(2),xi=xi';end
        x=xi(:,1);
        y=xi(:,2);
        polyorder=0;
        xlabelstring='x';
    end
end

if nargin==2, % Two arguments, might be x,y vectors or matrix and polyorder
    datasize=size(xi);
    if isscalar(yi),   %  matrix and polyorder
        if datasize(1)<datasize(2),xi=xi';end
        if isvector(xi),
            x=1:length(xi);
            xlabelstring='index number';
            y=xi';
            polyorder=yi;
        else
            x=xi(:,1);
            y=xi(:,2);
            polyorder=yi;
            xlabelstring='x';
        end
    else %  x,y vectors
        x=xi;
        y=yi;
        polyorder=0;
        xlabelstring='x';
    end
end

if nargin>2, % Three arguments, must be x,y vectors and polyorder
    x=xi;
    y=yi;
    xlabelstring='x';
end

x=reshape(x,[1,length(x)]);
y=reshape(y,[1,length(y)]);

if length(x)==length(y),
    if polyorder>0, % Plot and fit the data
        % Compute the fit
        [coef,S]=polyfit(x,y,polyorder);
        StdDevs=sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';
        % Plot the data and the polynomial fit
        xx=linspace(min(x),max(x)); % 100 values of x over range of x values
        yhat=polyval(coef,xx); % Predicted values of y for each xx
        plot(x,y,'.r',xx,yhat,'-b')
        if nargin==4,
            plot(x,y,datastyle,xx,yhat,'-b')
        end
        if nargin==5,
            plot(x,y,datastyle,xx,yhat,fitstyle)
        end
        axis([min(x) max(x) min(y) max(y)]);
        xlabel(xlabelstring);ylabel('y')
        title(['Number of data points= ' num2str(length(x)) ] )
        % Compute the correlation coefficient and R-Squared
        if IsOctave,
            cc=corr(polyval(coef,x'),y');
            RSquared=cc.^2;
        else
            cc=corrcoef(polyval(coef,x),y);
            RSquared=cc(2).^2;
        end %   if IsOctave,
        % Label the graph with the fit information
        text(min(x),max(y)-.05.*(max(y)-min(y)),['   Polynomial Order of fit = ' num2str(polyorder)] );
        text(min(x),max(y)-.1.*(max(y)-min(y)),['   Fit coefficients = ' num2str(coef)] );
        text(min(x),max(y)-.15.*(max(y)-min(y)),['   Standard deviation = ' num2str(StdDevs)] );
        text(min(x),max(y)-.2.*(max(y)-min(y)),['   R-Squared = ' num2str(RSquared)] );
    else % Just plot the data, no fit.
        plot(x,y)
        axis([min(x) max(x) min(y) max(y)]);
        xlabel(xlabelstring);ylabel('y')
        title(['Number of data points= ' num2str(length(x)) ] )
        coef=[];
        RSquared=[];
    end
else
    disp('Error: x and y must be the same size.')
    sizex=size(x)
    sizey=size(y)
end
    
if nargout==4, % Compute the bootsrtap
    BootstrapResultsMatrix=zeros(polyorder+1,NumTrials);
    clear bx by
    for trial=1:NumTrials,
        n=1;
        bx=x;
        by=y;
        while n<length(x)-1,
            if rand>.5,
                bx(n)=x(n+1);
                by(n)=y(n+1);
            end % if rand>.5,
            n=n+1;
        end % while n<length(xx)-1,
        coef=polyfit(bx,by,polyorder);
        BootstrapResultsMatrix(:,trial)=coef;
    end % for trial=1:NumTrials,
    disp(' ')
    disp('Bootstrap Results');
    BootstrapMean=mean(real(BootstrapResultsMatrix(:,:)'));
    BootstrapSTD=std(BootstrapResultsMatrix(:,:)');
    BootstrapIQR=IQrange(BootstrapResultsMatrix(:,:)');
    PercentRSD=100.*BootstrapSTD./abs(BootstrapMean);
    PercentIQR=100.*BootstrapIQR./abs(BootstrapMean);
    BootstrapMean=BootstrapMean(1:polyorder+1);
    BootstrapSTD=BootstrapSTD(1:polyorder+1);
    BootstrapIQR=BootstrapIQR(1:polyorder+1);
    PercentRSD=PercentRSD(1:polyorder+1);
    PercentIQR=PercentIQR(1:polyorder+1);
        disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
        disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
        disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
        disp(['Percent RSD:    ', num2str(PercentRSD)])
        disp(['Percent IQR:    ', num2str(PercentIQR)])
    BootResults(1,:)=BootstrapMean;
    BootResults(2,:)= BootstrapSTD ;
    BootResults(3,:)=PercentRSD;
    BootResults(4,:)=BootstrapIQR;
    BootResults(5,:)=PercentIQR;
end % if NumArgOut==3,

function b=IQrange(a)
% b = IQrange(a) returns the interquartile range of the values in a.  For
%  vector input, b is the difference between the 75th and 25th percentiles
%  of a.  For matrix input, b is a row vector containing the interquartile
%  range of each column of a. (c) T. C, O'Haver, 2012
% Example:
%  a=randn(10000,5);
%  IQrange(a)
% Divide by 1.34896 to get an estimate of the standard deviation withput
% outliers
mina=min(a);
sizea=size(a);
NumCols=sizea(2);
b=zeros(size(a));
for n=1:NumCols,
    b(:,n)=a(:,n)-mina(n);
end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));

function isOctave = IsOctave()
% Returns true if this code is being executed by Octave.
% Returns false if this code is being executed by MATLAB, or any other MATLAB
% variant.
% 
%    usage: isOctave = IsOctave()
    persistent octaveVersionIsBuiltIn;
    if (isempty(octaveVersionIsBuiltIn))
        octaveVersionIsBuiltIn = (exist('OCTAVE_VERSION', 'builtin') == 5);
        % exist returns 5 to indicate a built-in function.
    end
    isOctave = octaveVersionIsBuiltIn;
    % If OCTAVE_VERSION is a built-in function, then we must be in Octave.
    % Since the result cannot change between function calls, it is cached in a
    % persistent variable.  isOctave cannot be a persistent variable, because it
    % is the return value of the function, so instead the persistent result must
    % be cached in a separate variable.
%-----------------------------------------------------------------
function Processed=ProcessSignal(x,y,DerivativeMode,w,type,ends,Sharpen,factor1,factor2,SlewRate,MedianWidth)
% Command line function that performs smoothing and differentiation on the
% time-series data set x,y. Returns the processed signal as a vector that
% has the same shape as x, regardless of the shape of y.
% Version 1.0, April, 2013.  Tom O'Haver (toh@umd.edu)
%
% 'DerivativeMode' determines the derivative order (O, 1, 2, 3, 4, 5); 
%    See http://terpconnect.umd.edu/~toh/spectrum/Differentiation.html
%
% 'w' is the smooth width; 
% 'type' determines the smooth mode:
%        If type=0, the signal is not smoothed.
%        If type=1, rectangular (sliding-average or boxcar) 
%        If type=2, triangular (2 passes of sliding-average)
%        If type=3, pseudo-Gaussian (3 passes of sliding-average)
%        If type=4, Savitzky-Golay smooth 
% 'ends' controls how the "ends" of the signal (the first w/2 points and the last w/2 points) are handled.
%        If ends=0, the ends are zeroed
%        If ends=1, the ends are smoothed with progressively 
%        smaller smooths the closer to the end.
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html
%
% 'Sharpen' controls peak sharpening (resolution enhancement). 0=off; 1=on.
%   The sharpening str5ength is controled by the factor1 and factor2 The
%   optimum values depend on the peak shape and width. For details, see
%   http://terpconnect.umd.edu/~toh/spectrum/InteractiveResEnhance.htm). 
%  
% 'SlewRate' sets the maximum slew rate (maximum rate of change of the
%    signal when it is non-zero (useful for reducing the effect of large 
%    steps in the background).
% 
% 'MedianWidth' implements a median filter when it is non-zero, which
%    replaces each point in the signal with the median (rather than the
%    average) of MedianWidth adjacent points, cuseful for dealing with narrow
%    spike artifacts (set MedianWidth to the number of points in the spikes)
%
% EXAMPLE 1:
% >> ProcessSignal([1 3 5 7 9],[0 0 1 0 0],0,3,1,0,0,0,0,0,0)
% ans =
%             0      0.33333      0.33333      0.33333            0
% EXAMPLE 2:
% Use of smoothed second derivative to reduce influence of background.
% on the weak peak at x=5.
% >> x=[0:.01:10]';
% >> y=10./(1+(x/4).^2)+exp(-(x-5).^2)+.1.*randn(size(x));
% >> subplot(2,1,1);plot(x,y);
% >> subplot(2,1,2);plot(x,ProcessSignal(x,y,2,60,3,0,0,0,0,0,0))
%
% EXAMPLE 3: Single peak with random spikes. Use of median filter
% (MedianWidth=1) to remove single-point spikes.
% >> x=-5:.01:5;
% >> y=exp(-(x).^2);for n=1:1000,if randn()>2,y(n)=rand()+y(n);,end,end;
% >> plot(x,ProcessSignal(x,y,0,0,0,0,0,0,0,0,1))
%
if SlewRate,
    for n=1:length(y)-1,
        if y(n+1)-y(n)>SlewRate,y(n+1)=y(n)+SlewRate;end
        if y(n+1)-y(n)<-SlewRate,y(n+1)=y(n)-SlewRate;end
    end
end % SlewRate
if MedianWidth,
    mY=y;
    for n=1:length(x)-(1+MedianWidth*2),
        mY(n+MedianWidth)=median(y((n):(n+1+MedianWidth*2)));
        y=mY;
    end
end  % MedianWidth
if type==0,w=1;end
if type==4,
    if w<2,w=3;end
    if DerivativeMode>4,
        if w<5,w=5;end
    end 
    % The polynomial order, 2+DerivativeMode, must be less than the
    % frame size, 2*w+1, and 2*w+1 must be odd.  
        Processed=savitzkyGolayFilt(y,2+DerivativeMode,DerivativeMode,2*w+1);
        if DerivativeMode==1,Processed=-Processed;end
        if DerivativeMode==3,Processed=-Processed;end;
else
    switch DerivativeMode
        case 0
            Processed=fastsmooth(y,w,type,ends);
        case 1
            Processed=fastsmooth(deriv(x,y),w,type,ends);
        case 2
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            Processed=fastsmooth(D2,w,type,ends);
        case 3
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D3=fastsmooth(deriv(x,D2),w,1,ends);
            Processed=fastsmooth(D3,w,type,ends);
        case 4
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D4=fastsmooth(secderiv(x,D2),w,1,ends);
            D4=fastsmooth(D4,w,1,ends);
            
            Processed=fastsmooth(D4,w,type,ends);
        case 5
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D4=fastsmooth(secderiv(x,D2),w,1,ends);
            D4=fastsmooth(D4,w,1,ends);
            Processed=fastsmooth(deriv(x,D4),w,type,ends);
    end
end
if Sharpen,
    type=4; 
    if w<3;w=3;end
    Processed=enhance(x,Processed,factor1,factor2,w,type);
end
Processed=reshape(Processed,size(x));
% ----------------------------------------------------------------------
function Enhancedsignal=enhance(x,signal,factor1,factor2,SmoothWidth,type)
% Resolution enhancement function by even derivative method. The
% arguments factor1 and factor 2 are 2nd and 4th derivative weighting
% factors. Larger values of factor1 and factor2 will reduce the 
% peak width but will cause artifacts in the baseline near 
% the peak.  Adjust the factors for the the best compromise. 
% Use minimum smooth width needed to reduce excess noise. 
datasize=size(signal);
if datasize(1)>datasize(2),signal=signal';end
if type==4,
    Enhancedsignal = signal-factor1.*savitzkyGolayFilt(signal,4,2,2*SmoothWidth+1)+...
        factor2.*savitzkyGolayFilt(signal,6,4,2*SmoothWidth+1);
else
    d2=secderiv(x,signal);  % Computes second derivative
    d4=secderiv(x,d2);   % Computes fourth derivative
    Enhancedsignal = signal-factor1.*fastsmooth(d2,SmoothWidth,type)+...
        factor2.*fastsmooth(fastsmooth(fastsmooth(d4,SmoothWidth,2),SmoothWidth,2),SmoothWidth,2);
end
Enhancedsignal=Enhancedsignal';
% ----------------------------------------------------------------------
function d=secderiv(x,a)
% Second derivative of y with respect to x using 3-point central difference.
%  T. C. O'Haver, 2011.
n=length(a);
d=zeros(size(a));
for j = 2:n-2;
  x1=x(j-1);x2=x(j);x3=x(j+1);
  d(j)=((a(j+1)-a(j))./(x3-x2) - (a(j)-a(j-1))./(x2-x1))./((x3-x1)/2);
end
d(1)=d(2);
d(n)=d(n-1);
% ----------------------------------------------------------------------
function d=deriv(x,y)
% First derivative of y with respect to x using 2-point central difference.
%  T. C. O'Haver, 2011.
n=length(y);
d=zeros(size(y));
d(1)=(y(2)-y(1))./(x(2)-x(1));
d(n)=(y(n)-y(n-1))./(x(n)-x(n-1));
for j = 2:n-1;
  d(j)=(y(j+1)-y(j-1)) ./ (1.*(x(j+1)-x(j-1)));
end
% ----------------------------------------------------------------------
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
    case 0
       SmoothY=sa(Y,w,ends);  
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
  end% ----------------------------------------------------------------------
function y=savitzkyGolayFilt(x,N,DN,F,W,DIM)
%savitzkyGolayFilt Savitzky-Golay Filtering.
%   savitzkyGolayFilt(X,N,DN,F) filters the signal X using a Savitzky-Golay 
%   (polynomial) filter.  The polynomial order, N, must be less than the
%   frame size, F, and F must be odd.  DN specifies the differentiation
%   order (DN=0 is smoothing). For a DN higher than zero, you'll have to
%   scale the output by 1/T^DN to acquire the DNth smoothed derivative of
%   input X, where T is the sampling interval. The length of the input X
%   must be >= F.  If X is a matrix, the filtering is done on the columns
%   of X.
%
%   Note that if the polynomial order N equals F-1, no smoothing
%   will occur.
%
%   savitzkyGolayFilt(X,N,DN,F,W) specifies a weighting vector W with
%   length F containing real, positive valued weights employed during the
%   least-squares minimization. If not specified, or if specified as
%   empty, W defaults to an identity matrix.
%
%   savitzkyGolayFilt(X,N,DN,F,[],DIM) or savitzkyGolayFilt(X,N,DN,F,W,DIM)
%   operates along the dimension DIM.
%   Copyright (c) 2011, Diederick
%   See also savitzkyGolay, FILTER, sgolayfilt

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.4 $  $Date: 2009/08/11 15:47:54 $

error(nargchk(4,6,nargin,'struct'));

% Check if the input arguments are valid
if round(F) ~= F, error(generatemsgid('MustBeInteger'),'Frame length must be an integer.'), end
if rem(F,2) ~= 1, error(generatemsgid('SignalErr'),'Frame length must be odd.'), end
if round(N) ~= N, error(generatemsgid('MustBeInteger'),'Polynomial order must be an integer.'), end
if N > F-1, error(generatemsgid('InvalidRange'),'The Polynomial order must be less than the frame length.'), end
if DN > N, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

if nargin < 5 || isempty(W)
   % No weighting matrix, make W an identity
   W = ones(F,1);
else
   % Check for right length of W
   if length(W) ~= F, error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
   % Check to see if all elements are positive
   if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
end

if nargin < 6, DIM = []; end

% Compute the projection matrix B
pp = fix(-F./2):fix(F./2);
B = savitzkyGolay(pp,N,DN,pp,W);

if ~isempty(DIM) && DIM > ndims(x)
	error(generatemsgid('InvalidDimensions'),'Dimension specified exceeds the dimensions of X.')
end

% Reshape X into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

if size(x,1) < F, error(generatemsgid('InvalidDimensions'),'The length of the input must be >= frame length.'), end

% Preallocate output
y = zeros(size(x));

% Compute the transient on (note, this is different than in sgolayfilt,
% they had an optimization leaving out some transposes that is only valid
% for DN==0)
y(1:(F+1)/2-1,:) = fliplr(B(:,(F-1)/2+2:end)).'*flipud(x(1:F,:));

% Compute the steady state output
ytemp = filter(B(:,(F-1)./2+1),1,x);
y((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);

% Compute the transient off
y(end-(F+1)/2+2:end,:) = fliplr(B(:,1:(F-1)/2)).'*flipud(x(end-(F-1):end,:));

% Convert Y to the original shape of X
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end
% ----------------------------------------------------------------------
function [fc, df] = savitzkyGolay(x,n,dn,x0,W,flag)
% Function:
%       Savitzky-Golay Smoothing and Differentiation Filter
%       Copyright (c) 2011, Diederick
%       The Savitzky-Golay smoothing/differentiation filter (i.e., the
%       polynomial smoothing/differentiation filter, or  the least-squares
%       smoothing/differentiation filters) optimally fit a set of data
%       points to polynomials of different degrees. 
%       See for details in Matlab Documents (help sgolay). The sgolay
%       function in Matlab can deal with only symmetrical and uniformly
%       spaced data of even number.
%       This function presented here is a general implement of the sgolay
%       function in Matlab. The Savitzky-Golay filter coefficients for even
%       number, nonsymmetrical and nonuniformly spaced data can be
%       obtained. And the filter coefficients for the initial point or the
%       end point can be obtained too. In addition, either numerical
%       results or symbolical results can be obtained. Lastly, this
%       function is faster than MATLAB's sgolay.
%
% Usage:
%       [fc,df] = savitzkyGolay(x,n,dn,x0,flag)
%   input:
%       x    = the original data point, e.g., -5:5 
%       n    = polynomial order
%       dn   = differentation order (0=smoothing),  default=0
%       x0   = estimation point, can be a vector    default=0
%       W    = weight vector, can be empty          
%              must have same length as x0          default=identity
%       flag = numerical(0) or symbolical(1),       default=0
%
%   output:
%       fc   = filter coefficients obtained (B output of sgolay).
%       df   = differentiation filters (G output of sgolay).
% Notes:
% 1.    x can be arbitrary, e.g., odd number or even number, symmetrical or
%       nonsymmetrical, uniformly spaced or nonuniformly spaced, etc.       
% 2.    x0 can be arbitrary, e.g., the initial point, the end point, etc.
% 3.    Either numerical results or symbolical results can be obtained.
% Example:
%       sgsdf([-3:3],2,0,0,[],0)
%       sgsdf([-3:3],2,0,0,[],1)
%       sgsdf([-3:3],2,0,-3,[],1)
%       sgsdf([-3:3],2,1,2,[],1)
%       sgsdf([-2:3],2,1,1/2,[],1)
%       sgsdf([-5:2:5],2,1,0,[],1)     
%       sgsdf([-1:1 2:2:8],2,0,0,[],1)
% Author:
%       Diederick C. Niehorster <dcniehorster@hku.hk> 2011-02-05
%       Department of Psychology, The University of Hong Kong
%
%       Originally based on
%       http://www.mathworks.in/matlabcentral/fileexchange/4038-savitzky-golay-smoothing-and-differentiation-filter
%       Allthough I have replaced almost all the code (partially based on
%       the comments on the FEX submission), increasing its compatibility
%       with MATLABs sgolay (now supports a weight matrix), its numerical
%       stability and it speed. Now, the help is pretty much all that
%       remains.
%       Jianwen Luo <luojw@bme.tsinghua.edu.cn, luojw@ieee.org> 2003-10-05
%       Department of Biomedical Engineering, Department of Electrical Engineering
%       Tsinghua University, Beijing 100084, P. R. China  
% Reference
%[1]A. Savitzky and M. J. E. Golay, "Smoothing and Differentiation of Data
%   by Simplified Least Squares Procedures," Analytical Chemistry, vol. 36,
%   pp. 1627-1639, 1964.
%[2]J. Steinier, Y. Termonia, and J. Deltour, "Comments on Smoothing and
%   Differentiation of Data by Simplified Least Square Procedures,"
%   Analytical Chemistry, vol. 44, pp. 1906-1909, 1972.
%[3]H. H. Madden, "Comments on Savitzky-Golay Convolution Method for
%   Least-Squares Fit Smoothing and Differentiation of Digital Data,"
%   Analytical Chemistry, vol. 50, pp. 1383-1386, 1978.
%[4]R. A. Leach, C. A. Carter, and J. M. Harris, "Least-Squares Polynomial
%   Filters for Initial Point and Slope Estimation," Analytical Chemistry,
%   vol. 56, pp. 2304-2307, 1984.
%[5]P. A. Baedecker, "Comments on Least-Square Polynomial Filters for
%   Initial Point and Slope Estimation," Analytical Chemistry, vol. 57, pp.
%   1477-1479, 1985.
%[6]P. A. Gorry, "General Least-Squares Smoothing and Differentiation by
%   the Convolution (Savitzky-Golay) Method," Analytical Chemistry, vol.
%   62, pp. 570-573, 1990.
%[7]Luo J W, Ying K, He P, Bai J. Properties of Savitzky-Golay Digital
%   Differentiators, Digital Signal Processing, 2005, 15(2): 122-136.
%
%See also:
%       sgolay, savitzkyGolayFilt

% Check if the input arguments are valid and apply defaults
error(nargchk(2,6,nargin,'struct'));

if round(n) ~= n, error(generatemsgid('MustBeInteger'),'Polynomial order (n) must be an integer.'), end
if round(dn) ~= dn, error(generatemsgid('MustBeInteger'),'Differentiation order (dn) must be an integer.'), end
if n > length(x)-1, error(generatemsgid('InvalidRange'),'The Polynomial Order must be less than the frame length.'), end
if dn > n, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

% set defaults if needed
if nargin<6
    flag=false;
end
if nargin < 5 || isempty(W)
   % No weighting matrix, make W an identity
   W = eye(length(x0));
else
   % Check W is real.
   if ~isreal(W), error(generatemsgid('NotReal'),'The weight vector must be real.'),end
   % Check for right length of W
   if length(W) ~= length(x0), error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
   % Check to see if all elements are positive
   if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
   % Diagonalize the vector to form the weighting matrix
   W = diag(W);
end
if nargin<4
    x0=0;
end
if nargin<3
    dn=0;
end

% prepare for symbolic output
if flag
    x=sym(x);
    x0=sym(x0);
end

Nx  = length(x);
x=x(:);
Nx0 = length(x0);
x0=x0(:);

if flag
    A=ones(length(x),1);
    for k=1:n
        A=[A x.^k];
    end
    df = inv(A'*A)*A';                          % backslash operator doesn't work as expected with symbolic inputs, but the "slowness and inaccuracy" of this method doesn't matter when doing the symbolic version
else
    df = cumprod([ones(Nx,1) x*ones(1,n)],2) \ eye(Nx);
end
df = df.';

hx = [(zeros(Nx0,dn)) ones(Nx0,1)*prod(1:dn)];  % order=0:dn-1,& dn,respectively
for k=1:n-dn                                    % order=dn+1:n=dn+k
    hx = [hx x0.^k*prod(dn+k:-1:k+1)];
end

% filter coeffs
fc = df*hx'*W;