function [SineParams,my_std_peaks] = CBG_sineFit(x,y)
%Purpose: Estimation of noisy sine curve parameters by FFT and non linear fitting.
%
% Syntax:
%       [SineParams]=sineFit(x,y)
%       Input: x and y values, y=offs+amp+sin(2*pi*f*x+phi)+noise
%       Output: SineParams(1): offset (offs)
%               SineParams(2): amplitude (amp)
%               SineParams(3): frequency (f)
%               SineParams(4): phaseshift (phi)
%               my_std_peaks:  Standard deviation of the amplitude
%               difference of the peaks we find in the FFT. A high std
%               implies reliability.
%       yOut=offs+amp*sin(2*pi*f*x+phi)
%
% Example:
% % generate y(x)
% x=-4:5;
% y=1+2*(sin(2*pi*0.1*x+2)+0.3*randn(size(x)));%Sine + noise
% [SineP]=sineFit(x,y)
% figure;
% xx=x(1):(x(end)-x(1))/222:x(end);%better resolution
% plot(x,y,xx,SineP(1)+SineP(2)*sin(2*pi*SineP(3)*xx+SineP(4)));
% %uncomment following lines if you want to save y=f(x) and run it sineFitDemo
% %paramsClean=[1,2,0.1,2];
% % save('xy.mat','x','y','paramsClean');
%
%You may want to comment/uncomment the last statement (PlotResults) in the first function.
%Author: Peter Seibold
%% FFT
pi2=2*pi;
NumSamples=length(x);
T=x(2)-x(1);
fNy=1/(2*T);%Nyquist frequency
offs=mean(y);%DC value
y_m=y-offs;%FFT much better without offset
n = 128*2^nextpow2(NumSamples);%heavy zero padding
Y = fft(y_m,n);%Y(f)
n2=floor(n/2);
P2 = abs(Y/NumSamples);
P1 = P2(1:n2+1);
P1(2:end-1) = 2*P1(2:end-1);
fs = (0:n2)/n/T;% frequency scale
% %FFT parameters at peak
[maxFFT,maxFFTindx]=max(P1);%Peak magnitude and location
fPeak=fs(maxFFTindx);% f at peak
Phip=angle(Y(maxFFTindx))+pi/2;%Phi-Peak is for cos, sin(90°+alpha)=cos(betta), alpha=-betta
Phip=Phip-x(1)*fPeak*pi2;%shift for phi at x=0
%Better estimate for offset:
omega=pi2*fPeak;
offs=offs-maxFFT*(cos(omega*x(1)+Phip)-cos(omega*x(end)+Phip))/(omega*(x(end)-x(1)));
% Assess correctness of peaks
N_pks = 5;  % the maximum number of peaks we allow to find
[val_pks,idx_pks] = findpeaks(P1);
[val_pks_sort,tmp_idx] = sort(val_pks,'descend');
idx_pks_sort = idx_pks(tmp_idx);
plot(fs(idx_pks_sort(1:N_pks)),P1(idx_pks_sort(1:N_pks)),'linestyle','none','marker','*')
my_std_peaks = std(val_pks_sort(1:N_pks));
%% Fitting
paramsFFTp=[offs,maxFFT,fPeak,Phip];
if maxFFTindx<0.99*n2
  %FFT peak not at f-Nyquist
  NumPeaks=1;
  paramsFFT=paramsFFTp;
else
  %Samples per period close to 2, max FFT peak close to f-Nyquist
  %Set 1st evaluation point a little below f-Nyquist
  fIndxExtra1=round(maxFFTindx*.995);
  fExtra1=fs(fIndxExtra1);
  PhiExtra1=angle(Y(fIndxExtra1))+pi/2-x(1)*fExtra1*pi2;
  %extra f for evaluation left of max peak 
  fIndxExtra2=round(0.75*maxFFTindx);
  fExtra2=fs(fIndxExtra2);
  PhiExtra2=angle(Y(fIndxExtra2))+pi/2-x(1)*fExtra2*pi2;
  paramsFFT=[[offs,maxFFT,fPeak*.995,PhiExtra1];...
    [offs,0.8*maxFFT,fExtra2,PhiExtra2]]; 
  NumPeaks=2;
end
paramsOut=zeros(NumPeaks,6);%for regression outputs
%% find best fit in time domain
modelfun = @(paramc,x) paramc(1) + paramc(2) * sin(pi2*paramc(3)*x+paramc(4));
opts = statset('nlinfit');opts.MaxIter=1000;%620 is the limit in evaluated test set.
warning('off','all');%disable warnings from nlinfit
for i=1:NumPeaks
  [SineParams,~,~,~,MSE] = nlinfit(x,y,modelfun,paramsFFT(i,:),opts);
  %make frequency positive
  if SineParams(3)<0
    SineParams(3)=-SineParams(3);
    SineParams(4)=pi-SineParams(4);%sin(2*pi*-f-phi)=sin(2*pi*f+phi+pi)
  end
  %make amplitude positive
  if SineParams(2)<0
    SineParams(2)=-SineParams(2);
    SineParams(4)=SineParams(4)+pi;
  end
  paramsOut(i,:)=[SineParams,MSE,MSE];
  if NumSamples<5% && SineParams(3)<=fNy
    %No valid MSE from nlinfit if num samples <5
    %Overwrite MSE, set priority to 1st result (by *i)
    %will be overwritten again with max allowed amplitude
    paramsOut(i,5)=0.003*i;
  end
  if SineParams(3)>fNy
    %f larger than nyquist limit
    paramsOut(i,5)=Inf;%set MSE to terrible
  end
end
warning('on','all');

%% take best manipulated score
[MSEmin,MSEminIndx]=min(paramsOut(:,5));
SineParams=paramsOut(MSEminIndx,1:4);
%  Determine max allowed amplitude by MSEmin
if MSEmin<=0.00001 || ...%extremly good MSE
    NumSamples<5 || ... %no MSE with nlinfit for less than 5 samples
    (NumSamples==5 && SineParams(3)<0.8*paramsFFT(1,3)) ||... %num period propably <1
    (MSEmin<1 && x(end)-x(1)<0.5/SineParams(3))%propably less than 0.5 periods
  maxAmp=66*maxFFT;%max allowed amplitude
elseif MSEmin>0.3
  maxAmp=4*maxFFT;
elseif MSEmin>0.01
  maxAmp=6*maxFFT;
elseif MSEmin>0.001
  maxAmp=18*maxFFT;
else
  %very good MSE, 0.00001 < MSE <0.001
  maxAmp=33*maxFFT;
end
% maxAmp=0;%TEST! Force FFT output
if SineParams(2)>maxAmp || SineParams(3)>fNy
  %Best regression has too big amplitude or is over Nyquist limit,
  %take original FFT result
  SineParams=paramsFFTp;
  MSE=NaN;%for PlotResults
else
  MSE=paramsOut(MSEminIndx,6);%for PlotResults
end
%make phase between 0 and 2 pi
SineParams(4)=rem(SineParams(4),pi2);
if SineParams(4)<0
  SineParams(4)=SineParams(4)+pi2;
end

%Plot, uncomment following line or delete all following lines:
PlotResults(x,y,SineParams,paramsFFT,fs,P1,maxFFTindx,maxFFT,MSE);

%% Plot results (optional, uncomment statement above)
function PlotResults(x,y,SineParams,paramsFFT,fs,P1,maxFFTindx,maxFFT,MSE)
xstart=x(1);
xend=x(end);
x3b=(1:numel(x));
x3=(xend-xstart)/(numel(x)-1)*(x3b-1)+xstart;
x4b=1:0.01:numel(x);
x4=(xend-xstart)/(numel(x)-1)*(x4b-1)+xstart;
y5=SineParams(1)+SineParams(2)*sin(2*pi*SineParams(3)*x4+SineParams(4));%result
yFFT=paramsFFT(1,1)+paramsFFT(1,2)*sin(2*pi*paramsFFT(1,3)*x4+paramsFFT(1,4));
figure;%for time
plot(x3,y,'k.');%time series as dots
xlabel('Time [s]');
hold on;
pIn=plot(x3,y,'r-');%time series as line
pFFT=plot(x4,yFFT,'color',[0.9 0.9 0.9]);
pResult=plot(x4,y5,'b-');%result
legend([pIn,pResult,pFFT],'Input','Result', 'FFT peak');
hold off;
grid on;

figure;%for FFT
% title('FFT');
pFFTin=plot(fs,P1,'r-');
xlabel('Frequency [Hz]');
ylabel('Amplitude')
hold on;
pFFTmax=plot(fs(maxFFTindx),maxFFT,'r+','MarkerSize',12);%max FFT
pFFTresult=plot(SineParams(3),SineParams(2),'b+','LineWidth',2);
plot([SineParams(3),SineParams(3)],[0,max(max(P1)*1.01,SineParams(2))],'b-');
hLeg=legend([pFFTin,pFFTresult,pFFTmax],'Input',...
  ['Result:     ' num2str(SineParams(2),3) ', ' num2str(SineParams(3),3) ' Hz'],...
  ['max FFT:  ' num2str(maxFFT,3) ', ' num2str(fs(maxFFTindx),3) ' Hz'],...
  'Location','best');
title(hLeg,'        amplitude | frequency','FontSize',8);
hold off;
grid on;
disp(['Result:        y= ' num2str(SineParams(1)) ' + ' num2str(SineParams(2)) ...
  ' * sin(2*pi*' num2str(SineParams(3)) '+' num2str(SineParams(4)) ')   MSE: ' num2str(MSE)]);
disp(['FFT:           y= ' num2str(paramsFFT(1,1)) ' + ' num2str(paramsFFT(1,2)) ...
  ' * sin(2*pi*' num2str(paramsFFT(1,3)) '+' num2str(paramsFFT(1,4)) ')' ]);
