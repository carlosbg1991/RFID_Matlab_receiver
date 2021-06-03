function  [coef, RSquared,BootResults]=logplotfit(xi,yi,polyorder)
% logplotfit is a loglog version of plotfit; it accepts data in the form of
% a single vector, a pair of vectors, "x" and "y", or a 2xn or nx2 matrix
% with x in first row or column. Plots the log10(x) vs log10(y). If the
% input argument "polyorder" is supplied, it fits the log data to a
% polynomial of order "polyorder", plots the log data in red dots and the
% fit as a blue line, and displays the fit coefficients and R-squared in
% the upper left corner of the graph. Polyorder=1 for straight line, =2 for
% quadratic (parabola) etc. If the 3rd output argument (BootResults) is
% supplied, computes coefficient error estimates by the bootstrap method
% and returns the results in the matrix "BootResults" (of size 5 x
% polyorder+1). You can change the number of bootstrap samples in line 138. 
% Version 1, March 2014
% Examples:
% x=[.01 .1 1 10 100];y=[.02 .195 1.9 18.5 180];
% logplotfit(y); % plot y only vs index number, no fit.
% logplotfit(x,y); % plot x vs y only, data in separate vectors
% logplotfit([x;y]); % plot data only, data in matrix with x in first row
% logplotfit([x;y]'); % plot data only, data in matrix with x in first column
% logplotfit(y,2); % plot y vs index and fit to second order polynomial
% logplotfit(x,y,3); % plot x vs y and fit to third order polynomial
% logplotfit([x;y],3);  % plot and fit data in matrix
% [coef, RSquared]=logplotfit([x;y],2) % return fit coefficients and r-squared
% [coef, RSquared,BootResults]=logplotfit(x,y,polyorder) % computes
% coefficient error estimates by the bootstrap method and returns the
% results in the matrix "BootResults" 

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
% 
% process input arguments
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
            xlabelstring='log(x)';
        end
    else %  x,y vectors
        x=xi;
        y=yi;
        polyorder=0;
        xlabelstring='log(x)';
    end
end

if nargin==3, % Tree arguments, must be x,y vectors and polyorder
    x=xi;
    y=yi;
    xlabelstring='log(x)';
end

x=reshape(x,[1,length(x)]);
y=reshape(y,[1,length(y)]);

x=log10(x);
y=log10(y);

if length(x)==length(y),
    if polyorder>0, % Plot and fit the data
        % Compute the fit
        coef=polyfit(x,y,polyorder);
        
        % Plot the data and the polynomial fit
        xx=linspace(min(x),max(x));
        yhat=polyval(coef,xx);
        plot(x,y,'.r',xx,yhat,'-b')
        axis([min(x) max(x) min(y) max(y)]);
        xlabel(xlabelstring);ylabel('log(y)')
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
        text(min(x),max(y)-.15.*(max(y)-min(y)),['   R-Squared = ' num2str(RSquared)] );
    else % Just plot the data, no fit.
        plot(x,y)
        axis([min(x) max(x) min(y) max(y)]);
        xlabel(xlabelstring);ylabel('log(y)')
        title(['Number of data points= ' num2str(length(x)) ] )
        coef=[];
        RSquared=[];
    end
else
    disp('Error: x and y must be the same size.')
    sizex=size(x)
    sizey=size(y)
end
    
if nargout==3, % Compute the bootsrtap
    NumTrials=300; %<<< CHANGE number of bootstrap samples here >>>
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