% Demo of CLS regression of 5 random noisy overlapping Gaussian peaks
format compact 
clear
clf
noise=.01;
t=[1:5:1000];
amp=[rand(1) rand(1) rand(1) rand(1) rand(1)];  % Amplitudes of the peaks
pos=[200 300 400 600 800];   % Positions of the peaks
wid=[200 200 200 200 200];   % Widths of the peaks
% A = matrix containing one unit-amplidude peak in each of its rows
A = zeros(length(pos),length(t));
for k=1:length(pos)
  A(k,:)=lorentzian(t,pos(k),wid(k)); % You can use any peak function here
end
spectrum=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
  clf
  % Simulation
  ObservedSpectrum = spectrum + noise.*randn(size(t)); % adds noise
  hold on
  % Regress ObservedSpectrum onto matrix of reference spectra
  % Using straight multilinear least-squares regression formulation
  MeasuredAmp = ObservedSpectrum*A'*inv(A*A')
  % Alternative formulations
  % Using pseudoinverse
  % MeasuredAmp = ObservedSpectrum*pinv(A))';
  % Using matrix division
  % MeasuredAmp = ObservedSpectrum/A;  
  figure(1);plot(t,ObservedSpectrum,'m.')
    plot(t,amp(1).*A(1,:),t,amp(2).*A(2,:),...
      t,amp(3).*A(3,:),t,amp(4).*A(4,:),t,amp(5).*A(5,:))
  plot(t,MeasuredAmp(1).*A(1,:),t,MeasuredAmp(2).*A(2,:),...
      t,MeasuredAmp(3).*A(3,:),t,MeasuredAmp(4).*A(4,:),t,MeasuredAmp(5).*A(5,:))
    plot(t,MeasuredAmp(1).*A(1,:)+MeasuredAmp(2).*A(2,:)+MeasuredAmp(3).*A(3,:)+MeasuredAmp(4).*A(4,:)+MeasuredAmp(5).*A(5,:),'k')
  title('Analysis of mixture of peaks by multiple regression')
  title('Analysis of mixture of peaks by multiple regression')
  xlabel('Dots: Measured data; Colored lines true and measured component spectra')
  hold off
  % PercentErrors(k,:)=100.*(MeasuredAmplitudes-amp)
'    Peak 1    Peak 2    Peak 3    Peak 4    Peak 5'
TrueAmplitudes = amp
MeasuredAmplitudes = MeasuredAmp
Error = MeasuredAmplitudes-TrueAmplitudes;
RelativePercentError = 100.*Error./TrueAmplitudes