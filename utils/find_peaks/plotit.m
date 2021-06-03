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
% Tom O'Haver, Version 6.1; Version 6, March 2015, added optional specification of
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
% [coef, RSquared,StdDevs,BootResults]=plotit(x,y,polyorder) % computes
% coefficient error estimates by the bootstrap method and returns the
% results in the matrix "BootResults". Example of bootstrap calculation:
% >>x=1:100;;y=x+.5.*randn(size(x));
% >>[coef,RSquared,StdDevs,BootResults]=plotit(x,y,1);
% Bootstrap Results
% Mean:      0.99619     0.085319
% STD:       0.0013306    0.079314
% STD (IQR): 0.0013099    0.076073
% RSD:       0.133565      92.9623
% RSD (IQR): 0.131493       89.164
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

% Version 2, Copyright (c) 2017, Thomas C. O'Haver

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
NumTrials=100;  %   <<<< CHANGE if you wish
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
        title(['plotit.m     Number of data points= ' num2str(length(x)) ] )
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
        if polyorder==1;
            text(min(x),max(y)-.1.*(max(y)-min(y)),['   Slope = ' num2str(coef(1)) '   Intercept = ' num2str(coef(2))   ] );
        else
            text(min(x),max(y)-.1.*(max(y)-min(y)),['   Fit coefficients = ' num2str(coef)] );
        end
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
    disp('Plotit Error: x and y must be the same size.')
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
        disp(['Mean:        ', num2str(BootstrapMean)])
        disp(['STD:         ', num2str(BootstrapSTD)])
        disp(['STD (IQR):   ', num2str(BootstrapIQR)])
        disp(['% RSD:       ', num2str(PercentRSD)])
        disp(['% RSD (IQR): ', num2str(PercentIQR)])
    BootResults(1,:)=BootstrapMean;
    BootResults(2,:)= BootstrapSTD ;
    BootResults(3,:)=PercentRSD;
    BootResults(4,:)=BootstrapIQR;
    BootResults(5,:)=PercentIQR;
end % if NumArgOut==3,

function b=IQrange(a)
% b = IQrange(a) returns the estimated standard deviation (STD) computed
% from the interquartile range, divided by 1.34896, of the values in a. For
% vector input, the interquartile range is the difference between the 75th
%  and 25th percentiles of a. For matrix input, b is a row vector
%  containing the estimated STD of each column of a. (c) T. C, O'Haver,
%  2017
% Example:
%  a=randn(10000,5);
%  IQrange(a)

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
b=b/1.34896;

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