function [Height,Position,Width,BootResults]=bootgaussfit(x,y,plots)
% Fits a parabola (quadratic) to x vs ln(y), calculates the position,
% width, and height of the Gaussian from the three coefficients of the
% quadratic fit. Accurate only if the data have no baseline offset (that
% is, trends to zero far off the peak) and if there are no zeros or
% negative values in y. If plots=1, plots the raw data as red dots and the
% best-fit Gaussian as a line. If the 4th output argument (BootResults) is
% supplied, computes peak parameter error estimates and returns the results
% in the 3x5 matrix BootResults. You can change the number of bootstrap
% samples in line 39. 
% (c) Tom O'Haver, October 2012
%
% Example 1:
% x=50:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
%  [Height,Position,Width]=bootgaussfit(x,y,1) plots data and fit and
%   returns [Height,Position,Width] clustered around 100,100,100.
%
% Example 2:
% [Height,Position,Width,BootResults]=bootgaussfit(x,y,1);
%   As above but also computes peak parameter error estimates. If plots=0,
%   no output is generated but the table of error estimates is returned in
%   the matrix BootResults. If plots=1, data and fit are plotted and table
%   of error estimates is printed in command window.
%
% Example 3: As above, with 50% random noise, plots data set in blue and
% fitted region in red; compares result to peakfit results in figure 2.
%  x=-50:250;
%  y=100.*gaussian(x,100,100)+50.*randn(size(x));
%  figure(1);clg;plot(x,y,'.');
%  hold on;
%  [Height,Position,Width,BootResults]=bootgaussfit(x(100:200),y(100:200),1);
%    title('bootgaussfit.m results');hold off
%  figure(2);
%  [FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit([x;y]);
% 
ylim=max(y)/1000;
for p=1:length(y),
    if y(p)<ylim,y(p)=ylim;end % Prevents error stop if y=<0
end % for p=1:length(y),
z=log(y);
coef=polyfit(x,z,2);
a=coef(3);
b=coef(2);
c=coef(1);
Height=exp(a-c*(b/(2*c))^2);
Position=-b/(2*c);
Width=2.35482/(sqrt(2)*sqrt(-c));
if plots,
    xx=linspace(min(x), max(x));
    plot(x,y,'r.',xx,Height.*gaussian(xx,Position,Width))
    xlabel('x')
    ylabel('y')
end % if plots
if nargout==4, % Compute the bootsrtap
    NumTrials=300;
    if plots,disp('Computing bootstrap sampling statistics.....'),end
    BootstrapResultsMatrix=zeros(3,100);
    clear bx by
    bHeight=zeros(1,NumTrials);
    bPosition=zeros(1,NumTrials);
    bWidth=zeros(1,NumTrials);
    tic;
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
        [bHeight(trial),bPosition(trial),bWidth(trial)]=gaussfit(bx,by);
        BootstrapResultsMatrix(:,trial)=[bHeight(trial),bPosition(trial),bWidth(trial)];
    end % for trial=1:100,
    if plots,toc;end
    if plots,
        disp(' ')
        disp('               Height     Position       Width ');
    end % if plots
    BootstrapMean=mean(real(BootstrapResultsMatrix(:,:)'));
    BootstrapSTD=std(BootstrapResultsMatrix(:,:)');
    BootstrapIQR=IQrange(BootstrapResultsMatrix(:,:)');
    PercentRSD=100.*BootstrapSTD./BootstrapMean;
    PercentIQR=100.*BootstrapIQR./BootstrapMean;
    BootstrapMean=BootstrapMean(1:3);
    BootstrapSTD=BootstrapSTD(1:3);
    BootstrapIQR=BootstrapIQR(1:3);
    PercentRSD=PercentRSD(1:3);
    PercentIQR=PercentIQR(1:3);
    if plots,
        disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
        disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
        disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
        disp(['Percent RSD:    ', num2str(PercentRSD)])
        disp(['Percent IQR:    ', num2str(PercentIQR)])
    end % if plots
    BootResults(1,:)=BootstrapMean;
    BootResults(2,:)= BootstrapSTD ;
    BootResults(3,:)=PercentRSD;
    BootResults(4,:)=BootstrapIQR;
    BootResults(5,:)=PercentIQR;
end % if NumArgOut==4,

function [Height, Position, Width]=gaussfit(x,y)
% Converts y-axis to a log scale, fits a parabola
% (quadratic) to the (x,ln(y)) data, then calculates
% the position, width, and height of the
% Gaussian from the three coefficients of the
% quadratic fit.  This works only if the data have
% no baseline offset (that is, trends to zero far off the
% peak) and if there are no zeros or negative values in y.
% Example: [Height, Position, Width]=gaussfit([1 2 3],[1 2 1]) 
% returns Height = 2, Position = 2, Width = 2
maxy=max(y);
for p=1:length(y),
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
z=log(y);
coef=polyfit(x,z,2);
a=coef(3);
b=coef(2);
c=coef(1);
Height=exp(a-c*(b/(2*c))^2);
Position=-b/(2*c);
Width=2.35482/(sqrt(2)*sqrt(-c));

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);

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

