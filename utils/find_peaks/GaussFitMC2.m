function [MeasuredSTD ResultsMatrix]=GaussFitMC2(TrueHeight,TruePosition,TrueWidth,NumPoints,Noise,NumTrials)
% Monte Carlo simulation of the measurement of the peak height, position,
% and width of a noisy x,y Gaussian peak, comparing the gaussfit parabolic 
% fit to the fitgaussian iterative fit (included as subfunctions). 
% Plots the peak height, position, and width of the gaussfit results 
% vs the height, position, and width of the fitgaussian results, 
% plus a randomly selected single data set and Gaussian fit.
% The input arguments are:
% TrueHeight = True peak height of the Gaussian peak.
% TruePosition = True x-axis value at the peak maximum.
% TrueWidth =  True half-width (FWHM) of the Gaussian peak.
% NumPoints = Number of points taken for the least-squares fit
% Noise = standard deviaiton of (normally-distributed) random noise
% NumTrials = Number of trials used to compute standard deviations
% Optional Outputs:
% MeasuredSTD = computed standard deviation of Measured [gaussfitHeight
% gaussfitPosition gaussfitWidth fitgaussianHeight fitgaussianPosition fitgaussianWidth]
% ResultsMatrix = Measured [gaussfitHeight gaussfitPosition gaussfitWidth fitgaussianHeight fitgaussianPosition fitgaussianWidth] of each trial
% Example: GaussFitMC2(100,100,100,50,1,1000)
% ans =
%     0.2135    0.1403    0.4756    0.2103    0.1299    0.4590
global PEAKHEIGHTS
startx=TruePosition-TrueWidth/2;
endx=TruePosition+TrueWidth/2;
deltax=(endx-startx)/(NumPoints-1);
x=startx:deltax:endx;
ResultsMatrix=zeros(6,NumTrials);
for trial=1:NumTrials,
    NoiseArray=Noise*randn(size(x));
    y=TrueHeight*gaussian([startx:deltax:endx],TruePosition,TrueWidth)+NoiseArray;
    [Height, Position, Width]=gaussfit(x,y);
    start=[Position Width];    % generic initial guesses for peak position and width
    options = optimset('TolX',0.1);  % Determines how close the model must fit the data
    parameter=fminsearch(@(lambda)(fitgaussian(lambda,x,y)),start);
    ResultsMatrix(:,trial)=[Height, Position, Width, PEAKHEIGHTS, parameter(1), parameter(2)];
end
FullDataSetSTD=std(ResultsMatrix');
% Display results
MeasuredSTD=FullDataSetSTD;
subplot(2,2,1)
plotit(ResultsMatrix(1,:),ResultsMatrix(4,:),1);
title('Peak Heights');
xlabel('Parabolic fit')
ylabel('Interative fit')
axis([98.5 101.5 98.5 101.5])
subplot(2,2,2)
plotit(ResultsMatrix(2,:),ResultsMatrix(5,:),1);
title('Peak Positions');
xlabel('Parabolic fit')
ylabel('Interative fit')
axis([98.5 101.5 98.5 101.5])
subplot(2,2,3)
plotit(ResultsMatrix(3,:),ResultsMatrix(6,:),1);
title('Peak Widths');
xlabel('Parabolic fit')
ylabel('Interative fit')
axis([98.5 101.5 98.5 101.5])
subplot(2,2,4)
plot(x,y,'o',linspace(min(x),max(x)),Height*gaussian(linspace(min(x),max(x)),Position,Width))
title('Single data set');
xlabel('x')
ylabel('y')

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
z=log(y);
coef=polyfit(x,z,2);
a=coef(3);
b=coef(2);
c=coef(1);
Height=exp(a-c*(b/(2*c))^2);
Position=-b/(2*c);
Width=2.35482/(sqrt(2)*sqrt(-c));

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