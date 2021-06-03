% Demonstration of effect of random noise on the power transform method by
% taking the 1/n power of the measured heights and areas of data raised to
% the power n. Effect of n on precision and goodness of fit. Takes several
% seconds to complete.
clear P
n(power)=power;
Noise=.1; % STD of random white noise
height=1; % Original peak height when power=1
width=4; % Original peak width when power=1
x=-5:.05:5;
for power=1:9 % Power
    for trial=1:100
        y=(height.*gaussian(x,0,width)+Noise.* randn(size(x))).^ n(power);
        [results,GOF]=peakfit([x;y],0,0,0,0,0,0,0,0,0,0);
        P(trial,:)=results;
        E(trial,:)=GOF;
        a(trial)=trapz(x,y);
        ap(trial)=trapz(x,y).^(1/n(power));
        w(trial)=halfwidth(x,y);
    end
    MeanHeight(n(power))=mean(P(:,3));
    RSDHeight(n(power))=rsd(P(:,3));
    Error(n(power))=mean(E(:,1));
    R2(n(power))=mean(E(:,2));
    CorrectedHeight(n(power))= mean(P(:,3)).^(1/n(power));
    CorrectedArea(n(power))= mean(P(:,5)).^(1/n(power));
end
figure(1)
clg
semilogy(n,MeanHeight,n,RSDHeight,n,Error,n,R2,n,CorrectedHeight,n,CorrectedArea,'LineWidth',2)
xlabel('power')
title('Blue=MeanHeight Green=RSDHeight  Red=%FitError   Cyan=R2   Magenta=CorrectedHeight  Yellow=CorrectedArea ')
grid