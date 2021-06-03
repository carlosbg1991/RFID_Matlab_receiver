% Demonstration of the limits of Classical Least Squares (multilinear
% regression) of two very closely-spaced "noiseless" overlapping Gaussian
% peaks caused by the numerical precision of the computer software. This
% demonstrates three different mathematical formulation of the
% least-squares calculation that give different results when the numerical
% precision limits of the computer are reached. But practically, the
% difference between these methods is unlikey to be seen; even the tiniest
% bit of added random noise (line 15) or signal instability produces a
% far greater error.
format compact
format long g
clear
clf
ADCbits=12;
noise=0.0; % Also try small anounts of noise like .00000001;
t=[1:10:1000];
amp=[1 .9];  % Amplitudes of the peaks
pos=[500 511];   % Positions of the peaks
wid=[300 300];   % Widths of the peaks
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
  % Regress ObservedSpectrum onto matrix of reference spectra. Can be done
  % bt several different mathematical methods that give the same results.
  % Using straight multilinear least-squares regression formulation
  tic
  MeasuredAmpInverse = ObservedSpectrum*A'*inv(A*A');
InverseTime=toc 
  % Using pseudoinverse
  tic
  MeasuredAmpPseudoinverse = ObservedSpectrum*pinv(A);
  PseudoinverseTime=toc
  % Using matrix division
  tic
  MeasuredAmpMatrixDivision = ObservedSpectrum/A;  
  MatrixDivisionTime=toc
  figure(1);plot(t,ObservedSpectrum,'m.')
    plot(t,amp(1).*A(1,:),t,amp(2).*A(2,:))
  plot(t,MeasuredAmpInverse(1).*A(1,:),t,MeasuredAmpInverse(2).*A(2,:))
    plot(t,MeasuredAmpInverse(1).*A(1,:)+MeasuredAmpInverse(2).*A(2,:),'k')
  title('Analysis of mixture of peaks by multiple regression')
  xlabel('Dots: Measured data; Colored lines are true and measured component spectra (blue and green)')
  hold off
disp(' ')
SeparationToWidthRatio=(pos(2)-pos(1))./mean(wid)
disp(['Method 1, the inverse method:  MeasuredAmp=ObservedSpectrum*A''*inv(A*A'')'] )
disp(MeasuredAmpInverse)
ErrorInverse = MeasuredAmpInverse-amp;
RelativePercentErrorInverse = 100.*ErrorInverse./amp
disp('--------------------------------------------------------------------')

disp('Method 2, the pseudo-inverse method:  MeasuredAmp=ObservedSpectrum*pinv(A)' )
disp(MeasuredAmpPseudoinverse)
ErrorPseudoinverse = MeasuredAmpPseudoinverse-amp;
RelativePercentErrorPseudoinverse = 100.*ErrorPseudoinverse./amp
disp('--------------------------------------------------------------------')

disp('Method 2, the matrix division method:  MeasuredAmp=ObservedSpectrum/A' )
disp(MeasuredAmpMatrixDivision)
ErrorMatrixDivision = MeasuredAmpMatrixDivision-amp;
RelativePercentErrorMatrixDivision = 100.*ErrorMatrixDivision./amp