function BootstrapIterativeFit(TrueHeight,TruePosition,TrueWidth,NumPoints,Noise,NumTrials)
% Bootstrap estimation of the standard deviation and interquartile range of
% an iterative least-squares fit to a single noisy Gaussian peak.  Plots
% the histograms of peak height, position, and width of the bootstrap samples.
%
% The input arguments are:
% TrueHeight = True peak height of the Gaussian peak.
% TruePosition = True x-axis value at the peak maximum.
% TrueWidth =  True half-width (FWHM) of the Gaussian peak.
% NumPoints = Number of points taken for the least-squares fit
% Noise = standard deviation of (normally-distributed) random noise
% NumTrials = Number of bootstrap samples
%
% Example:
% >> BootstrapIterativeFit(100,100,100,20,10,100);
%        Peak Height  Peak Position  Peak Width
% mean:   99.27028    100.4002       94.5059
% STD:    2.8292      1.3264         2.9939
% IQR:    4.0897      1.6822         4.0164
% IQR/STD Ratio:    1.3518
%
global PEAKHEIGHTS

% Create x,y data set with noise
startx=TruePosition-TrueWidth;
endx=TruePosition+TrueWidth;
deltax=(endx-startx)/(NumPoints-1);
x=startx:deltax:endx;
NoiseArray=Noise*randn(size(x));
y=TrueHeight*gaussian(x,TruePosition,TrueWidth)+NoiseArray;
% y=fastsmooth(y,5,1,1); % Optional smoothing

% Perform bootstrap
BootstrapResultsMatrix=zeros(3,NumTrials);
clear xx yy
cutoff=0.5;
for trial=1:NumTrials,
    n=1;
    xx=x;
    yy=y;
    while n<length(x)-1,
        if rand>cutoff,
            xx(n)=x(n+1);
            yy(n)=y(n+1);
        end
        n=n+1;
    end
    start=[TruePosition TrueWidth];    % generic initial guesses for peak position and width
    % options = optimset('TolX',0.1);  % Determines how close the model must fit the data
    % Modified form to work in both Matlab abd Octave:
    parameter=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),start);
    ThisResult=[PEAKHEIGHTS parameter(1) parameter(2)];
    BootstrapResultsMatrix(:,trial)=ThisResult;
end
BootstrapMean=mean(real(BootstrapResultsMatrix'));
BootstrapSTD=std(BootstrapResultsMatrix');
BootstrapIQR=IQrange(BootstrapResultsMatrix');
BootstrapMedian=median(BootstrapResultsMatrix');
BootstrapMode=mode(BootstrapResultsMatrix');

% Display results
disp('       Peak Height  Peak Position  Peak Width')
disp(['mean:   ' num2str(BootstrapMean)]);
disp(['STD:    ' num2str(BootstrapSTD)]);
disp(['IQR:    ' num2str(BootstrapIQR)]);
disp(['IQR/STD Ratio:    ' num2str(mean(BootstrapIQR./BootstrapSTD))]);
disp(['Median:    ' num2str(BootstrapMedian)]);
disp(['Mode:    ' num2str(BootstrapMode)]);
% Plot original data set and histograms
subplot(2,2,1)
plot(x,y,'o',linspace(min(x),max(x)),BootstrapMean(1)*gaussian(linspace(min(x),max(x)),BootstrapMean(2),BootstrapMean(3)))
title('Original data set');
xlabel('x')
ylabel('y')

subplot(2,2,2)
hist(real(BootstrapResultsMatrix(1,:)),10);
title('Peak Height');

subplot(2,2,3)
hist(real(BootstrapResultsMatrix(2,:)),10);
title('Peak Position');

subplot(2,2,4)
hist(real(BootstrapResultsMatrix(3,:)),10);
title('Peak Width');

function err = fitgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal.
global PEAKHEIGHTS
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);

function b=IQrange(a)
% b = IQrange(a) returns the interquartile range divided by 1.34896, of the
% values in a, which is an estimate of the standard deviation of a without
% ouliers. If a is a vector, b is the difference between the 75th and 25th
% percentiles of a, divided by 1.34896. If b is a matrix, b is a row vector
% containing the interquartile range of each column of a, divided by
% 1.34896. The factor 1.34896 makes this quantity equal on average to the
% standard deviation (SD) of normal distributions, and thus IQrange is a
% better estimate of the standard deviation without ouliers for a normal
% distribution.
%  T. C. O'Haver, 2017
% Example:
%  a=randn(10000,5);
%  IQrange(a)
% Divide by 1.34896 to get an estimate of the standard deviation withput
% outliers
mina=min(a);
sizea=size(a);
NumCols=sizea(2);
b=zeros(size(a));
for n=1:NumCols
    b(:,n)=a(:,n)-mina(n);
end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));
b=b./1.34896;