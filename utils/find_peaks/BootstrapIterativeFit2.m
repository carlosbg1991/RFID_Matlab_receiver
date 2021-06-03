function BootstrapIterativeFit2(TrueHeight1,TruePosition1,TrueWidth1,TrueHeight2,TruePosition2,TrueWidth2,NumPoints,Noise,NumTrials)
% Bootstrap estimation of the standard deviation and interquartile range
% of an iterative least-squares fit to two overlapping noisy Gaussian
% peaks. Plots the histograms of peak height, position, and width of the 
% bootstrap samples.
%
% The input arguments are:
% TrueHeight1 = True peak height of Gaussian peak 1.
% TruePosition1 = True x-axis value at the peak maximum 1.
% TrueWidth1 =  True half-width (FWHM) of Gaussian peak 1.
% TrueHeight2 = True peak height of Gaussian peak 2.
% TruePosition2 = True x-axis value at peak maximum 2.
% TrueWidth2 =  True half-width (FWHM) of Gaussian peak 2.
% NumPoints = Number of points taken for the least-squares fit
% Noise = standard deviation of (normally-distributed) random noise
% NumTrials = Number of bootstrap samples
%
% Example:
% >> BootstrapIterativeFit2(100,100,100,50,200,100,100,10,200);
% Full Fit Results from NumTrials signals
% Peak 1    Height     Position     Width
% mean:   97.98754      100.1033      99.62307
% STD:    15.288       7.7692      13.8835
% IQR:    7.85432      7.17797      12.4954
% Peak 2    Height     Position     Width
% mean:   51.16718       199.049      99.59502
% STD:    9.73134      19.0192      28.3057
% IQR:    8.64353      13.6295      21.4718
%  
% Bootstrap Results from one signal
% Peak 1    Height     Position     Width
% mean:   79.7318      95.9337      90.6087
% STD:    15.5817      7.54896      11.3722
% IQR:    15.9845       6.9505      11.1164
% Peak 2    Height     Position     Width
% mean:   55.59346      183.7168      126.9655
% STD:    10.5914      19.3139      19.7621
% IQR:    10.3563      19.4515      18.5433
%
global PEAKHEIGHTS

% Create x,y data set with noise
startx=TruePosition1-2*TrueWidth1;
endx=TruePosition2+2*TrueWidth2;
deltax=(endx-startx)/(NumPoints-1);
x=startx:deltax:endx;
FullFitResultsMatrix=zeros(6,NumTrials);
for trial=1:NumTrials,
    NoiseArray=Noise*randn(size(x));
    y=TrueHeight1*gaussian(x,TruePosition1,TrueWidth1)+TrueHeight2*gaussian(x,TruePosition2,TrueWidth2)+NoiseArray;
    start=[TruePosition1 TrueWidth1 TruePosition2 TrueWidth2];    % generic initial guesses for peak position and width
    options = optimset('TolX',0.1);  % Determines how close the model must fit the data
    parameter=fminsearch(@(lambda)(fitgaussian(lambda,x,y)),start,options);
    % fminsearch(@(lambda)(fitgaussian(lambda,x,y)),newstart,options);
    ThisResult=[PEAKHEIGHTS(1) parameter(1) parameter(2) PEAKHEIGHTS(2) parameter(3) parameter(4)];
    FullFitResultsMatrix(:,trial)=ThisResult;
end
FullFitMean=mean(real(FullFitResultsMatrix'));
FullFitSTD=std(FullFitResultsMatrix');
FullFitIQR=IQrange(FullFitResultsMatrix');
% Display results
figure(1)
xx=linspace(min(x),max(x));
plot(x,y,'.',x,FullFitMean(1)*gaussian(x,FullFitMean(2),FullFitMean(3)),xx,FullFitMean(4)*gaussian(xx,FullFitMean(5),FullFitMean(6)))
title('Original data set and fitted peaks');
xlabel('Green = Peak 1     Red = Peak 2 ')
ylabel('y')
disp(['Full Fit Results from ' num2str(NumTrials) ' signal measurements'])
disp('Peak 1    Height     Position     Width')
disp(['mean:   ' num2str(FullFitMean(1:3))]);
disp(['STD:    ' num2str(FullFitSTD(1:3))]);
disp(['IQR:    ' num2str(FullFitIQR(1:3))]);
disp(['IQR/STD Ratio:    ' num2str(mean(FullFitIQR(1:3)./FullFitSTD(1:3)))]);
disp('Peak 2    Height     Position     Width')
disp(['mean:   ' num2str(FullFitMean(4:6))]);
disp(['STD:    ' num2str(FullFitSTD(4:6))]);
disp(['IQR:    ' num2str(FullFitIQR(4:6))]);
disp(['IQR/STD Ratio:    ' num2str(mean(FullFitIQR(4:6)./FullFitSTD(4:6)))]);

% Perform bootstrap
BootstrapResultsMatrix=zeros(6,NumTrials);
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
    start=[TruePosition1 TrueWidth1 TruePosition2 TrueWidth2];    % generic initial guesses for peak position and width
    options = optimset('TolX',0.1);  % Determines how close the model must fit the data
    parameter=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),start,options);
   
    ThisResult=[PEAKHEIGHTS(1) parameter(1) parameter(2) PEAKHEIGHTS(2) parameter(3) parameter(4)];
    BootstrapResultsMatrix(:,trial)=ThisResult;
end
BootstrapMean=mean(real(BootstrapResultsMatrix'));
BootstrapSTD=sqrt(2)*IQrange(BootstrapResultsMatrix');
BootstrapIQR=IQrange(BootstrapResultsMatrix');
% Display results
disp(' ');
disp('Bootstrap Results from one signal measurement')
disp('Peak 1    Height     Position     Width')
disp(['mean:   ' num2str(BootstrapMean(1:3))]);
disp(['STD:    ' num2str(BootstrapSTD(1:3))]);
disp(['IQR:    ' num2str(BootstrapIQR(1:3))]);
disp(['IQR/STD Ratio:    ' num2str(mean(BootstrapIQR./BootstrapSTD))]);
disp('Peak 2    Height     Position     Width')
disp(['mean:   ' num2str(BootstrapMean(4:6))]);
disp(['STD:    ' num2str(BootstrapSTD(4:6))]);
disp(['IQR:    ' num2str(BootstrapIQR(4:6))]);
disp(['IQR/STD Ratio:    ' num2str(mean(BootstrapIQR./BootstrapSTD))]);

% Plot original data set and histograms


figure(2)
subplot(3,2,1)
hist(real(BootstrapResultsMatrix(1,:)),10);
title('Peak Height 1');

subplot(3,2,2)
hist(real(BootstrapResultsMatrix(2,:)),10);
title('Peak Position 1');

subplot(3,2,3)
hist(real(BootstrapResultsMatrix(3,:)),10);
title('Peak Width 1');

subplot(3,2,4)
hist(real(BootstrapResultsMatrix(4,:)),10);
title('Peak Height 2');

subplot(3,2,5)
hist(real(BootstrapResultsMatrix(5,:)),10);
title('Peak Position 2');

subplot(3,2,6)
hist(real(BootstrapResultsMatrix(6,:)),10);
title('Peak Width 2');


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