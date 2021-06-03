% Demo of "classical" CLS regression .
% Example of multiple unknown samples; there are 10 samples, each
% containing 5 components, and the spectra have 100 wavelengths. In this
% case the 10 samples are repeat measurements of the same unknown mixture
% (but the same method can be used if they are different mixtures). Each
% time you run the simulation, you get a different mixture of the 5
% components.
format short g
clear
clf
t=[1:10:1000]; % Defines the "number of wavelengths"
noise=.01;
NumWavelengths=length(t);
NumSamples=10; % Defines the "number of samples"
% Each component is represented by a single Gaussian peak at a different
% wavelength, overlapping with adjuacent peaks.
% The number of components is defined by the length of amp, pos, and wid.
amp=[rand(1) rand(1) rand(1) rand(1) rand(1)];  % Amplitudes of the 5 peaks
% amp=[5 4 3 2 1];
pos=[200 300 400 600 800];   % Positions of the 5 peaks
wid=[200 200 200 200 200];   % Widths of the 5 peaks
% A = matrix containing one unit-amplitude peak in each of its rows
A = zeros(length(pos),length(t));
for k=1:length(pos)
    A(k,:)=gaussian(t,pos(k),wid(k)); % You can use any peak function here
end
spectrum=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
clf
% Simulation. NumSamples=number of unknown samples to be measured
for s=1:NumSamples,
    % In this case, the samples are repeat measurtememtns of one mixture
    ObservedSpectrum(s,:) = spectrum + noise.*randn(size(t));
end
hold on
% Regress ObservedSpectrum onto matrix of reference spectra. Can be done
% bt several different mathematical methods that give the same results.
% Using straight multilinear least-squares regression formulation
disp('   Component 1  Component 2  Component 3  Component 4  Component 5')
MeasuredAmp = ObservedSpectrum*A'*inv(A*A')
% Using pseudoinverse
% MeasuredAmp = ObservedSpectrum*pinv(A);
% Using matrix division
% MeasuredAmp = ObservedSpectrum/A;
figure(1);plot(t,ObservedSpectrum,'.')
plot(t,amp(1).*A(1,:),t,amp(2).*A(2,:),...
    t,amp(3).*A(3,:),t,amp(4).*A(4,:),t,amp(5).*A(5,:))
plot(t,MeasuredAmp(:,1)*A(1,:),t,MeasuredAmp(:,2)*A(2,:),...
    t,MeasuredAmp(:,3)*A(3,:),t,MeasuredAmp(:,4)*A(4,:),t,MeasuredAmp(:,5)*A(5,:))
%  plot(t,MeasuredAmp(:,1)*A(1,:)+MeasuredAmp(:,2)*A(2,:)+MeasuredAmp(:,3)*A(3,:)+MeasuredAmp(:,4)*A(4,:)+MeasuredAmp(:,5)*A(5,:),'k')
title('Analysis of multiple unknown mixtures of known components by multiple regression')
xlabel('Colored dots: Measured spectra of unknown mixtures; Colored lines: measured components')
hold off
% Compute the relative percent difference between the measured and true
% amplitudes
PercentErrors=100.*(mean(MeasuredAmp)-amp)./amp