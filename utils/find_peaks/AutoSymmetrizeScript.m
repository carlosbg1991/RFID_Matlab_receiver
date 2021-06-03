% function EMGplusfirstderivative
% Demonstration of symmetricalization of an exponentially modified
% Gaussian (EMG) by weighted addition of the first derivative, weighting
% factor s/(lambda*deltat^2). Shows that the peak width is reduced, the
% peak made more symmetrical, and the area is preserved.
format compact
format short g
clear
A=1;
s=50; % standard deviation of Gaussian peak
tau=100; % Exponential broadening factor
mu=1000; % peak center
deltat=5; % sampling interval
Noise=0.01; % random white noise
SmoothWidth=41; % Smooth width to reduce effect of noise
t=0:deltat:2500; % time axis
EMG=A.*expgaussian2(t,mu,s,tau)+Noise.*randn(size(t));
EMG=fastsmooth(EMG,SmoothWidth,2);

[hw,slopeafter1,slopeafter2]=halfwidth(t,EMG,mu);
r=slopeafter1/slopeafter2;
EstimatedTau=r.*-hw/4-hw/8;
mu=t(val2ind(EMG,max(EMG)));
FirstDeriv=derivxy(t,EMG);
% InitialWidth=halfwidth(t,EMG,mu);
trial=0;
StartFactor=EstimatedTau/2;
EndFactor=EstimatedTau*2;
for factor=StartFactor:1:EndFactor
   trial=trial+1;
   % Weighted sum of EMG and its first derivative
   % EMGSym = symmetricalized EMG
   EMGSym=EMG+factor.*FirstDeriv;
   EMGSym=fastsmooth(EMGSym,SmoothWidth,2);
   Miny(trial)=min(EMGSym);
   TrialTau(trial)=factor;
   MinD=min(FirstDeriv);
   EMGmin(trial)=EMGSym(val2ind(FirstDeriv,MinD));
   [~,slopeafter1,slopeafter2]=halfwidth(t,EMGSym,mu);
   SlopeRatioAfter(trial)=slopeafter1/slopeafter2;
   % disp([trial factor EMGmin(trial) SlopeRatioAfter(trial)])
end

coef=polyfit(SlopeRatioAfter,TrialTau,3);
EstimateBySymmetry=polyval(coef,-1);

for tt=length(TrialTau):-1:1
    if Miny(tt)>.1*min(Miny)
        EstimateByZeroBaseline=TrialTau(tt);
        break
    end
end

if exist('EstimateByZeroBaseline','var')   
else   
    EstimateByZeroBaseline=0;
end

figure(2)
clf
coef=plotit(SlopeRatioAfter,TrialTau,3);
grid
ylabel('TrialTau (time units)')
xlabel('Ratio of leading edge to trailing edge slope')
title(['Ratio equals -1 when peak is symmetrical, at tau = '  num2str(EstimateBySymmetry) '.'])

figure(3)
clf
plot(TrialTau,Miny)
xlabel('TrialTau (time units)')
 title(['Graph hits zero at tau = '  num2str(EstimateByZeroBaseline) '.'])
ylabel('Minimum of EMGSym')


figure(1);
clf
factor=EstimateBySymmetry;
EMGSym=EMG+factor.*FirstDeriv;
EMGSym=fastsmooth(EMGSym,SmoothWidth,2);
plot(t,EMG,t,EMGSym)
xlabel(['Blue: Original exponentially broadened peak.   Red: AutoSymmetrized peak.  Estimated tau: ' num2str(EstimateBySymmetry) '.'])
title('Automated determination of exponential broadening factor, tau')
disp(' ')

 disp('         Noise      SmoothWidth     tau   EstBySymmetry  EstByZeroBaseline ')
disp([Noise SmoothWidth tau EstimateBySymmetry EstimateByZeroBaseline])
% 
