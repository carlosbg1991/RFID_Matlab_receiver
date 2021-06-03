% Script that compares standard deviation of slope and intercept 
% for a first-order least-squares fit computed by Monte Carlo
% simulation of 1000 repeats to predictions made by closed-form
% algebraic equations. Change the values in lines 5-10 to suit.
% Version 2, April 2015, changed STD function to be compatible with Octave.
%
Noise=5;        % Standard deviation of the noise in the data
TrueSlope=1;     % Actual underlying slope of the data
TrueIntercept=0; % Actual underlying y-intercept of the data
minx=1;          % minimum value of x
maxx=100;         % maximum value of x
NumPoints=100;    % Number of points in the data
%
% Monte Carlo simulation
FullResults=zeros(10000,2);
for repeat=1:1000 
    x=minx:(maxx-minx)/NumPoints:maxx;
    NoiseArray=Noise*(randn(size(x))); % Constant noise added to each data point
    y=TrueIntercept+x*TrueSlope+NoiseArray;
    PredSet=polyfit(x,y,1);
    if repeat==1,
        FullResults=PredSet;
    else
        FullResults=[FullResults; PredSet];
    end
end % repeat
FullMean=mean(FullResults);
if IsOctave,
    FullSTD=[stdev(FullResults(:,1)) stdev(FullResults(:,2))];
else
    FullSTD=std(FullResults);
end
%
% Algebraic least-squares computation.
Sxx=sum((x-mean(x)).^2);
Syy=sum((y-mean(y)).^2);
Sxy=sum((x-mean(x)).*(y-mean(y)));
CalcSlope=Sxy./Sxx;
CalcIntercept=mean(y)-CalcSlope*mean(x);
Sy=sqrt((Syy-CalcSlope^2*Sxx)/(NumPoints-2));
SDSlope=Sy/sqrt(Sxx);
SDIntercept=Sy*sqrt(1./(NumPoints-(sum(x).^2)./sum(x.^2)));
PredMean=[CalcSlope CalcIntercept];
PredSTD=[SDSlope SDIntercept];
%
% Results comparison
disp(' ')
if IsOctave,
   disp(['NumPoints = ', num2str(NumPoints) '     SD Noise = ', num2str(stdev(NoiseArray)) '    x-range = ' num2str(maxx-minx)])
else
   disp(['NumPoints = ', num2str(NumPoints) '     SD Noise = ', num2str(std(NoiseArray)) '    x-range = ' num2str(maxx-minx)]) 
end
disp('                    MC Simulation           Algebraic Prediction')
disp('                    slope     intercept     slope  intercept')
disp(['Mean:               ', num2str([FullMean PredMean]) ])
disp(['Standard deviation: ', num2str([FullSTD PredSTD]) ])
disp('--------------------------------------------------------------')
plot(x,y,'o',linspace(min(x),max(x)),TrueSlope*linspace(min(x),max(x))+TrueIntercept)
xlabel('x')
ylabel('y')