function ShapeTestA(x,y,BaselineCorrectionMode)
% Fits data in x,y to 6 different asymmetrical function shapes and compares
% the fitting errors and R2 values. Plots each trial fit in a separate
% figure window. Requires peakfit.m in the path.
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
%   x=1:100;y=expgaussian(x,40,10,-10)';ShapeTestA(x,y,0)
%
% Example 2: Loads complex saved signal, tests one selected peak.
% (requires testsignalS.mat in the Matlab/Octave path)
% load testsignalA.mat
% ShapeTestA(x(820:1000),y(820:1000),1)
%
% Related function: ShapeTestS (for symmetrical peak shapes)
disp(' ')
disp('Peak Shape           fitting error   R2')
% 'One  Gaussian'
figure(1)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,1,0,10,0,BaselineCorrectionMode,0);
disp(['One Gaussian:           ' num2str(MeanFitError) ])
subplot(2,1,1)
title('One unconstrained Gaussian')

% 'Two Gaussians'
figure(2)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,2,1,0,10,0,BaselineCorrectionMode,0); 
disp(['Two Gaussians:          ' num2str(MeanFitError) ])
subplot(2,1,1)
title('Two unconstrained Gaussians')

% Exponentially modified Gaussian, variable time constant
figure(3)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,31,0,10,0,BaselineCorrectionMode,0);
disp(['Exp. modified Gaussian: ' num2str(MeanFitError) ])
subplot(2,1,1)
title('One exponentially broadened Gaussian, variable time constant')

% 'BiGaussian'
figure(4)
asymmetry=1.5;
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,14,asymmetry,10,0,BaselineCorrectionMode,0);
disp(['BiGaussian:             ' num2str(MeanFitError) ])
subplot(2,1,1)
title('One BiGaussian, fixed asymmetry')

% 'Exponential pulse'
figure(5)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,9,0,20,0,BaselineCorrectionMode,0);
disp(['Exponential pulse:      ' num2str(MeanFitError) ])
subplot(2,1,1)
title('One unconstrained exponential pulse')

% 'Alpha function'
figure(6)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,19,0,20,0,BaselineCorrectionMode,0);
disp(['Alpha function:         ' num2str(MeanFitError) ])
subplot(2,1,1)
title('One unconstrained alpha function')
