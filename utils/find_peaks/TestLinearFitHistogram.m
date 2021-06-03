% Script that compares standard deviation of slope and intercept for a
% first-order least-squares fit computed in three ways: (1) Monte Carlo
% simulation of NumTrials repeats line 14; (2)closed-form algebraic
% equations; and (3) simple bootstrap. Displays the histograms of the
% standard deviations of the repeat trials in figure 1 (slope) and in
% figure 2 (intercept). Prints table of mean and standard deviation of the
% repeat trials. Make NumTrials and NumRepeats larger for greater precision
Noise=10;
TrueSlope=2;
TrueIntercept=10;
minx=20;
maxx=50;
NumPoints=25;
NumTrials=50; % Make larger for greater precision
NumRepeats=50; % Make larger for greater precision
for RepeatTrial=1:NumTrials,
    FullResults=zeros(10000,2);
    % Prediction using empirical equation derived from Monte Carlo simulations
    % disp(['3*log(pi)*Noise./((maxx-minx)*(sqrt(NumPoints))) = '
    % num2str(3*log(pi)*Noise./((maxx-minx)*(sqrt(NumPoints))))])
    for repeat=1:NumRepeats;
        x=minx:(maxx-minx)/NumPoints:maxx;
        NoiseArray=Noise*(randn(size(x)));
        y=TrueIntercept+x*TrueSlope+NoiseArray;
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
    
    clear BootstrapResultsMatrix xx yy
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
        PolyFitResults=polyfit(xx,yy,1);
        BootstrapResultsMatrix(trial,:)=[PolyFitResults(1) PolyFitResults(2)]';
    end
    %subplot(2,2,1)
    BootstrapMean=mean(real(BootstrapResultsMatrix));
    BootstrapSTD=std(BootstrapResultsMatrix);
    if RepeatTrial==1,
       RepeatResults=[FullSTD PredSTD BootstrapSTD];
    else
       RepeatResults=[RepeatResults; [FullSTD PredSTD BootstrapSTD]];
    end
end % RepeatTrial

figure(1)
subplot(2,2,1)
hist(RepeatResults(:,1),10)
title('Monte Carlo simulation')
xlabel('Standard deviation of slope')
%axis([.05 .15 0 50])
subplot(2,2,2)
hist(RepeatResults(:,3),10)
axis([.05 .15 0 50])
title('Algebraic')
xlabel('Standard deviation of slope')
subplot(2,2,3)
hist(RepeatResults(:,5),10)
%axis([.05 .15 0 50])
title('Bootstrap')
xlabel('Standard deviation of slope')

figure(2)
subplot(2,2,1)
hist(RepeatResults(:,2),10)
title('Monte Carlo simulation')
xlabel('Standard deviation of Intercept')
%axis([2.5 5.5 0 50])
subplot(2,2,2)
hist(RepeatResults(:,4),10)
%axis([2.5 5.5 0 50])
title('Algebraic')
xlabel('Standard deviation of Intercept')
subplot(2,2,3)
hist(RepeatResults(:,6),10)
%axis([2.5 5.5 0 50])
title('Bootstrap')
xlabel('Standard deviation of Intercept')

MeanResults=mean(RepeatResults);
STDResults=std(RepeatResults);
disp(['NumPoints = ', num2str(NumPoints) '     SD Noise = ', num2str(std(NoiseArray)) '    x-range = ' num2str(maxx-minx)])
disp('    Simulation        Algebraic equation    Bootstrap method')
disp('    SDslope  SDint      SDslope  SDint      SDslope  SDint')
disp(MeanResults)
disp(STDResults);

