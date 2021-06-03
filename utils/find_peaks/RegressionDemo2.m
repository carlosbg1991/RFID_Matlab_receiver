% Demo of CLS regression of random noisy overlapping Gaussian peaks
% version 2 compares simple and weighted regression
format compact 
clear
clf
noise=.001;
t=[1:10:1000];
amp=[rand(1) rand(1) rand(1) rand(1) rand(1)];  % Amplitudes of the peaks
amp=amp./2;
pos=[200 300 400 600 800];   % Positions of the peaks
wid=[400 400 400 400 400];   % Widths of the peaks
% A = matrix containing one unit-amplidude peak in each of its rows
A = zeros(length(pos),length(t));
for k=1:length(pos)
  A(k,:)=gaussian(t,pos(k),wid(k)); % You can use any peak function here
end
spectrum=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
  clf
  % Simulation
  ObservedSpectrum = spectrum + noise.*randn(size(t)); % adds noise
  hold on
  % Regress ObservedSpectrum onto matrix of reference spectra
  % Using multilinear least-squares regression formulation
  % w=10.^-ObservedSpectrum;
  w=ObservedSpectrum;
  MeasuredAmp = ObservedSpectrum/A;  
  MeasuredAmp2=([w' w' w' w' w']' .* A)'\(ObservedSpectrum .* w)';
  MeasuredAmp2=MeasuredAmp2';
  % Using matrix inverse
  %   MeasuredAmp = ObservedSpectrum*A'*inv(A*A');
  % Using pseudoinverse
  %   MeasuredAmp = ObservedSpectrum*pinv(A))';
  % Using matrix division
  %   MeasuredAmp = ObservedSpectrum/A;  
  figure(1);plot(t,ObservedSpectrum,'m.')
    plot(t,amp(1).*A(1,:),t,amp(2).*A(2,:),...
      t,amp(3).*A(3,:),t,amp(4).*A(4,:),t,amp(5).*A(5,:))
  plot(t,MeasuredAmp(1).*A(1,:),t,MeasuredAmp(2).*A(2,:),...
      t,MeasuredAmp(3).*A(3,:),t,MeasuredAmp(4).*A(4,:),t,MeasuredAmp(5).*A(5,:))
    plot(t,MeasuredAmp(1).*A(1,:)+MeasuredAmp(2).*A(2,:)+MeasuredAmp(3).*A(3,:)+MeasuredAmp(4).*A(4,:)+MeasuredAmp(5).*A(5,:),'k')
  title('Analysis of mixture of peaks by multiple regression')
  xlabel('Dots: Measured data; Colored lines true and measured component spectra')
  ylabel('Measured absorbance')
  hold off
  % PercentErrors(k,:)=100.*(MeasuredAmplitudes-amp)
'    Peak 1    Peak 2    Peak 3    Peak 4    Peak 5'
TrueAmplitudes = amp
MeasuredAmplitudes = MeasuredAmp
MeasuredAmplitudes2 = MeasuredAmp2
Error = MeasuredAmplitudes-TrueAmplitudes;
Error2 = MeasuredAmplitudes2-TrueAmplitudes;
RelativePercentError = 100.*Error./TrueAmplitudes
RelativePercentError2 = 100.*Error2./TrueAmplitudes
MeanError=mean(abs(RelativePercentError))
MeanError2=mean(abs(RelativePercentError2))
