function ShapeTestS(x,y,BaselineCorrectionMode)
% Fits data in x,y to 7 different symmetrical function shapes and compares
% fitting errors abd R2 values. Plots each trial fit in a separate figure
% window. Requires peakfit.m in the path.
% BaselineCorrectionMode can be 0, 1, 2 or 3:
% 0 (default) does not subtract baseline from data segment; 
% 1 interpolates a linear baseline from the edges of the data 
%  segment and subtracts it from the signal (assumes that the 
%  peak returns to the baseline at the edges of the signal); 
% 2 is like mode 1 except that it computes a quadratic curved baseline; 
% 3 compensates for a flat baseline without reference to the 
%  signal itself (best if the peak does not return to the
%  baseline at the edges of the signal).
%
% Example 1: Single isolated peak, no background or noise.
%  x=1:100;y=voigt(x,50,5,1);ShapeTestS(x,y,0)
%
% Example 2: Loads complex saved signal, tests one selected peak.
% (requires testsignalS.mat in the Matlab/Octave path)
%  load testsignalS.mat
%  ShapeTestS(x(70:170),y(70:170),1)
%
% Related function: ShapeTestA (for asymmetrical peak shapes)
warning off all

% Single Gaussian
figure(1)
[FitResults1,MeanFitError1]=peakfit([x' y'],0,0,1,1,0,10,0,BaselineCorrectionMode,0);
subplot(2,1,1)
title('Single unconstrained Gaussian')

% Two Gaussians
figure(2)
[FitResults2,MeanFitError2]=peakfit([x' y'],0,0,2,1,0,10,0,BaselineCorrectionMode,0); 
subplot(2,1,1)
title('Two unconstrained Gaussians')

% Single Lorentzian
figure(3)
[FitResults3,MeanFitError3]=peakfit([x' y'],0,0,1,2,0,10,0,BaselineCorrectionMode,0);
subplot(2,1,1)
title('Single Lorentzian')

% logistic distribution
figure(4)
[FitResults4,MeanFitError4]=peakfit([x' y'],0,0,1,3,0,10,0,BaselineCorrectionMode,0);
subplot(2,1,1)
title('Logistic distribution')
% variable-alpha Voigt
figure(5)
[FitResults5,MeanFitError5]=peakfit([x' y'],0,0,1,30,0,10,0,BaselineCorrectionMode,0);

subplot(2,1,1)
title('variable-alpha Voigt profile')
% 'Variable shape Pearson'
figure(6)
[FitResults6,MeanFitError6]=peakfit([x' y'],0,0,1,32,0,10,0,BaselineCorrectionMode,0);
subplot(2,1,1) 
title('Variable shape Pearson')

% 'Variable Gaussian/Lorentzian blend
figure(7)
[FitResults7,MeanFitError7]=peakfit([x' y'],0,0,1,33,0,10,0,BaselineCorrectionMode,0);
subplot(2,1,1)
title('Variable Gaussian/Lorentzian blend')

disp(' ')
disp('Peak Shape            fitting error   R2')
disp(['Single Gaussian:          ' num2str(MeanFitError1) ])
disp(['Two Gaussians:            ' num2str(MeanFitError2) ])
disp(['Single Lorentzian:        ' num2str(MeanFitError3) ])
disp(['Logistic distribution:    ' num2str(MeanFitError4) ])
disp(['Voigt profile:            ' num2str(MeanFitError5) ])
disp(['Single Pearson:           ' num2str(MeanFitError6) ])
disp(['Gauss/Lorentz blend:      ' num2str(MeanFitError7) ])