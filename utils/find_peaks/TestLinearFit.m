% Script that compares standard deviation of slope and intercept for a
% first-order least-squares fit computed by random-number simulation to
% predictions made by closed-form algebraic equations and to a bootstrap.
NumPoints=100; % Number of data points in signal
Noise=10; % Standard deviation of the random noise
TrueSlope=2;
TrueIntercept=10;
minx=1;
maxx=30;
SmoothWidth=0;
FullResults=zeros(100,2);
% Prediction using empirical equation derived from Monte Carlo simulations
% disp(['3*log(pi)*Noise./((maxx-minx)*(sqrt(NumPoints))) = '
% num2str(3*log(pi)*Noise./((maxx-minx)*(sqrt(NumPoints))))])
disp(['NumPoints = ', num2str(NumPoints) '     SD Noise = ', num2str(Noise) '    x-range = ' num2str(maxx-minx)])
disp('    Simulation                  Algebraic equation        Bootstrap method')
disp('    SDslope        SDint        SDslope      SDint        SDslope       SDint')
for test=1:10,
    for repeat=1:100
        x=minx:(maxx-minx)/NumPoints:maxx;
        % Choose one of the noise models below
        NoiseArray=Noise*(randn(size(x))); % normal white constant noise
        % NoiseArray=Noise*bluenoise(x); % Blue noise
        % NoiseArray=Noise*pinknoise(length(x)); % Pink noise        
        % NoiseArray=bimodal(x,Noise,30,-30); % bimodal white constant noise
        % NoiseArray=(TrueIntercept+x*TrueSlope).*(randn(size(x)))./10; % normal white proportional noise
        y=TrueIntercept+x*TrueSlope+NoiseArray;
        if SmoothWidth,
            y=fastsmooth(y,SmoothWidth,3,1);
        end
        PredSet=polyfit(x,y,1);
        if repeat==1,
            FullResults=PredSet;
        else
            FullResults=[FullResults; PredSet];
        end
    end % repeat
    FullMean=mean(FullResults);
    FullSTD=std(FullResults);
    % Least-squares computation without polyfit function.
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
    BootstrapResultsMatrix=zeros(100,2);
    clear BootstrapResultsMatrix xx yy
    cutoff=0.5;
    for trial=1:500,
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
        PolyFitResults=polyfit(xx,yy,1);
        BootstrapResultsMatrix(trial,:)=[PolyFitResults(1) PolyFitResults(2)]';
    end
    %subplot(2,2,1)
    SubsetMean=mean(real(BootstrapResultsMatrix));
    SubsetSTD=sqrt(2)*std(BootstrapResultsMatrix);
    if SmoothWidth,disp(['Smooth Width = ', num2str(SmoothWidth)]);end
    disp([FullSTD PredSTD SubsetSTD]);
    figure(1)
    plot(x,y,'o',linspace(min(x),max(x)),TrueSlope*linspace(min(x),max(x))+TrueIntercept)
    figure(2)
    plot(xx,yy,'o',linspace(min(xx),max(xx)),TrueSlope*linspace(min(xx),max(xx))+TrueIntercept)
end
disp('--------------------------------------------------------------')