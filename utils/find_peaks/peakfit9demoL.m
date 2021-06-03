% Demonstration script for peakfit.m version 9 (Requires modelpeaks.m and
% peakfit.m in the Matlab path). Creates a noisy model signal of three
% peaks of known shapes, positions, and widths, but unknown heights.
% Compares multilinear regression (shape 50) in Figure 1 with unconstrained
% iterative non-linear least squares in Figure 2. For the multilinear
% regression mode (shape 50), the FIXEDPARAMETERS matrix must be supplied
% listing the peak shape (column 1), position (column 2), and width (column
% 3) of each peak, one row per peak. For example, in this demo
% FIXEDPARAMETERS =
%      2   400    70
%      2   500    85
%      2   560    90
%
% Although the fitting error is greater (and the R2 lower) for multilinear
% regression (shape 50), the errors in the measured peak heights are often
% less.
% T. C. O'Haver 2018
% Related script: peakfit9demo.m is similar but uses Gaussian peaks
% (specified in the FIXEDPARAMETERS matrix and in the PeakShape vector)
format short g
format compact
% User settable values
x=250:1:800;
PeakShape=[2 2 2]; % 1=Gaussian,  2=Lorentzian, etc
PeakHeights=[10 30 20];
PeakPositions=[400 500 560];
PeakWidths=[70 85 90];
extra=[0 0 0]; % Must match length of PeakShape even if values are zeros.
Noise=.5;
BaselineMode=0;

% Calculate simulated signal with noise
NumPeaks=length(PeakHeights);
model=modelpeaks(x,NumPeaks,PeakShape,PeakHeights,PeakPositions,PeakWidths,extra);
y=model+Noise.*randn(size(x)); % Add random noise

% Compute fit by multilinear regression
FIXEDPARAMETERS=[PeakShape(1) PeakPositions(1) PeakWidths(1);PeakShape(2) PeakPositions(2) PeakWidths(2);PeakShape(3)   PeakPositions(3)    PeakWidths(3)];
figure(1)
clf
[FitResults,GOF]=peakfit([x;y],0,0,3,50,0,1,0,BaselineMode,FIXEDPARAMETERS);
disp('Multilinear regression results (known position and width):')
disp('           Peak     Position      Height         Width     Area')
disp(FitResults)
MeanHeightError=100.*mean(abs((PeakHeights-FitResults(:,3)')./PeakHeights))';
disp(['% fitting error=' num2str(GOF(1)) '    R2= ' num2str(GOF(2))  '    % MeanHeightError= ' num2str(MeanHeightError)  ] )
disp(' ')
% Compute fit by iterative non-linear least squares
figure(2)
clf
NumTrials=5;
start=[PeakPositions(1) PeakWidths(1) PeakPositions(2) PeakWidths(2)  PeakPositions(3) PeakWidths(3)];
[FitResults,GOF]=peakfit([x;y],0,0,3,PeakShape,extra,NumTrials,start,BaselineMode);
MeanHeightError=100.*mean(abs((PeakHeights-FitResults(:,3)')./PeakHeights))';
disp('Unconstrained iterative non-linear least squares results:')
disp('           Peak     Position      Height         Width     Area')
disp(FitResults)
disp(['% fitting error=' num2str(GOF(1)) '    R2= ' num2str(GOF(2))  '    % MeanHeightError= ' num2str(MeanHeightError)  ] )


disp(' ')