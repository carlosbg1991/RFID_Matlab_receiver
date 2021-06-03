function [MeasuredSTD PredictedSTD ResultsMatrix]=GaussFitMC(TrueHeight,TruePosition,TrueWidth,NumPoints,Noise,NumTrials)
% Monte Carlo simulation of the measurement of the peak height, position,
% and width of a noisy x,y Gaussian peak, by fitting a parabola 
% to a region of x vs ln(y) equal to one Gaussian half-width, 
% using the gaussfit function (included as a subfunction). Plots a 20-bin
% histogran of the NumTrials results of peak height, position,
% and width, plus a randomly selected single data set and Gaussian fit.
% The input arguments are:
% TrueHeight = True peak height of the Gaussian peak.
% TruePosition = True x-axis value at the peak maximum.
% TrueWidth =  True half-width (FWHM) of the Gaussian peak.
% NumPoints = Number of points taken for the least-squares fit
% Noise = standard deviaiton of (normally-distributed) random noise
% NumTrials = Number of trials used to compute standard deviation
% Optional Outputs:
% MeasuredSTD = computed standard deviation of [Height Position Width]
% PredictedSTD = Predicted standard deviation (SD) based on empirical
% equations based on a similar Monte Carlo simulation: 
%    SD of the peak height = 1.73*noise/sqrt(NumPoints),
%    SD of the peak position = noise/sqrt(NumPoints,) 
%    SD of the peak width = 3.62*noise/sqrt(NumPoints). 
% ResultsMatrix = Measured [Height Position Width] of each trial
% Example: GaussFitMC(100,100,100,50,1,1000)
% ans =
%     0.2454    0.1365    0.5182
startx=TruePosition-TrueWidth/2;
endx=TruePosition+TrueWidth/2;
deltax=(endx-startx)/(NumPoints-1);
x=startx:deltax:endx;
ResultsMatrix=zeros(3,NumTrials);
for trial=1:NumTrials,
    NoiseArray=Noise*randn(size(x));
    y=TrueHeight*gaussian([startx:deltax:endx],TruePosition,TrueWidth)+NoiseArray;
    [Height, Position, Width]=gaussfit(x,y);
    ResultsMatrix(:,trial)=[Height, Position, Width];
end
FullDataSetSTD=std(ResultsMatrix');
RootFactor=sqrt(length(x));
HeightFactor=1.73;
PositionFactor=.99;
WidthFactor=3.62;
% Display results
TypicalHeightPositionWidth=[Height, Position, Width];
MeasuredSTD=FullDataSetSTD;
PredictedSTD=[Noise.*HeightFactor./RootFactor Noise.*PositionFactor./RootFactor Noise.*WidthFactor./RootFactor];
subplot(2,2,1)
hist(ResultsMatrix(1,:),20)
title('Histogram of Peak Heights');
xlabel('Height')
ylabel('Probability')
subplot(2,2,2)
hist(ResultsMatrix(2,:),20)
title('Histogram of Peak Positions');
xlabel('Position')
ylabel('Probability')
subplot(2,2,3)
hist(ResultsMatrix(3,:),20)
title('Histogram of Peak Widths');
xlabel('Width')
ylabel('Probability')
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

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);