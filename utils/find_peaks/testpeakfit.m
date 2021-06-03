% Test script for peakfit.m version 7.7 or later. Requires peakfit.m,
% ExpBroaden.m, gaussian.m, lorewtzian.m, humps.m, DataMatrix3.mat to be
% in the path. Version 6, December 2015
warning off all
disp('-------------------------------------------------------------------')
disp('Test script for peakfit.m version 7.7 or later. Requires peakfit.m,')
disp('ExpBroaden.m, gaussian.m, lorentzian.m, humps.m, DataMatrix3.mat')
disp('to be in the path. Make sure you have the latest peakfit.m; get it')
disp(' from http://tinyurl.com/cey8rwh ')
disp(' ')
echo on
% Example 1:  Signal is a matrix: Fits x vs y data with a single Gaussian
%  peak model. 
%
%       Peak number  Peak position   Height     Width      Peak area
x=[0:.1:10];y=exp(-(x-5).^2);[Results,FitError]=peakfit([x' y'])
echo off
drawnow;echo on
%---------------------------------------------------------------------
% Example 2: Signal is a single vector:  Create a small data set and fit to
%  a Gaussian model
[Results,FitError]=peakfit([0 0 0 0 0 0 0 0 0 0 1 2 4 6 7 6 4 2 1 0 0 0 0 0 0 0 0 0 0])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 3:   Like Example 1, except that random noise is added to the y
%  data. Compare output to Example 1.   
x=[0:.1:10];y=exp(-(x-5).^2)+.1*randn(1,length(x));[Results,FitError]=peakfit([x' y'])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%--------------------------------------------------------------------- 
% Example 4:   Fits a noisy two-peak signal with a double Gaussian model
% (NumPeaks=2).
%
x=[0:.1:10];y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(1,length(x));
[Results,FitError]=peakfit([x' y'],5,19,2,1,0,1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 5: Attempt to fit a portion of the humps function, 0.7 units wide
% and centered on x=0.3, with a single (NumPeaks=1) Pearson function
% (peakshape=4) with extra=3 (controls shape of Pearson function).
%
x=[0:.005:1];y=humps(x);[Results,FitError]=peakfit([x' y'],.3,.7,1,4,3)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
  %---------------------------------------------------------------------  
% Example 6: Creates a data matrix 'smatrix', fits a portion to a two-peak
%  Gaussian model, takes the best of 10 trials.  Returns optional output
%  arguments Results and FitError.
%
x=[0:.005:1];y=(humps(x)+humps(x-.13)).^3;smatrix=[x' y'];
[Results,FitError]=peakfit(smatrix,.4,.7,2,1,0,10)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 7: As above, but specifies the first-guess position and width of
%  the two peaks, in the order [position1 width1 position2 width2]
%
[Results,FitError]=peakfit([x' y'],.4,.7,2,1,0,10,[.3 .1 .5 .1])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 8: As above, returns the vector x1 containing 600 interploated
%  x-values for the model peaks and the matrix y1 containing the y values
%  of each model peak at each xi.  Type plot(xi,yi(1,:)) to plot peak 1 or
%  plot(xi,yi) to plot all peaks.
%
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit(smatrix,.4,.7,2,1,0,10);
clf
plot(xi,yi)
title('Plot of model peaks evaluated at 600 x-values')
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
% --------------------------------------------------------------------
% Example 9:  Like Example 7, but sets autozero=1 in the last argument.
%  Default is autozero=0.
%
[Results,FitError]=peakfit([x' y'],.4,.7,2,1,0,10,[.3 .1 .5 .1],1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 10: Fits a group of three peaks near x=2400 in DataMatrix3 with
% three Gaussians (shape=1). The residuals are very wavy, suggesting that
% model is not quite right.
%
load DataMatrix3 
[Results,FitError]=peakfit(DataMatrix3,2420,440,3,1,31,1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on

% Fitting the same data with three exponentially-broadened Gaussians
% (shape=5) takes longer but yields a much lower fit error and the residuls
% are a little more random and less wavy:
disp('Calculating........')
[Results,FitError]=peakfit(DataMatrix3,2420,440,3,5,31,1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 11:  Example of an unstable fit to a signal consisting of two
%  Gaussian peaks of equal height (1.0). The peaks are too highly
%  overlapped for a stable fit, even though the fit error is small and the
%  residuals are unstructured. Each time you re-generate this signal, it
%  gives a different fit, with the peaks heights varying about 15% from
%  signal to signal.  
%
x=[0:.1:10]';y=exp(-(x-5.5).^2)+exp(-(x-4.5).^2)+.01*randn(size(x));
[Results,FitError]=peakfit([x y],5,19,2,1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%
% More stable results can be obtained using the fixed-width Gaussian
% model, but that is justified only if the experiment is legitimately 
% expected to yield peaks of known width. Even though the fit error is
% slightly greater, the peak heights are actually more accurate.

[Results,FitError]=peakfit([x y],5,19,2,11,0,10,0,0,[1.66 1.66])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 12: Accuracy of autozero, for a single noiseless Gaussian on a
%  curved background, specifying center (5.0) and window (5.5), but using
%  placeholders (zeros) to accept the default values for NumPeaks,
%  peakshape, extra, NumTrials, and start (Version 2.5 and above). The
%  "true" values of peak parameters are 5.00, 1.00, 1.66, and 1.77.
%
%  No baseline correction (ATOZERO=0)
x=[0:.1:10]';y=1./(1+x.^2)+exp(-(x-5).^2);[Results,FitError]=peakfit([x y],5,5.5,0,0,0,0,0,0)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
% 
% Linear autozero (autozero=1)
x=[0:.1:10]';y=1./(1+x.^2)+exp(-(x-5).^2);[Results,FitError]=peakfit([x y],5,5.5,0,0,0,0,0,1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
% 
% Quadratic autozero (autozero=2): 
x=[0:.1:10]';y=1./(1+x.^2)+exp(-(x-5).^2);[Results,FitError]=peakfit([x y],5,5.5,0,0,0,0,0,2)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%
% Flat baseline correction (autozero=3)
x=[0:.1:10]';y=1./(1+x.^2)+exp(-(x-5).^2);[Results,FitError,Baseline]=peakfit([x y],5,5.5,0,0,0,0,0,3)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 13: Same as example 4, fit with fixed-width Gaussian (shape 11),
%  width=1.666. (added in version 2.6)
% 
x=[0:.1:10];y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
[Results,FitError]=peakfit([x' y'],0,0,2,11,0,0,0,0,[1.666 1.666]);
Results
FitError
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 14: Peak area measurements. Same as the example in Figure 15 on
%  Integration and Peak Area Measurement.  All four peaks have the same
%  theoretical peak area (1.772). The four peaks can be fit together in one
%  fitting operation using a 4-peak Gaussian model, with only rough
%  estimates of the first-guess positions and widths.  The peak areas thus
%  measured are much more accurate than the perpendicular drop method:
%
x=0:.05:18;
y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2);
[Results,FitError]=peakfit([x;y],0,0,4,1,0,1,[4 2 9 2 12 2 14 2],0,0)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
% 
% This works well even in the presence of substantial amounts of random
% noise:
%
x=0:.05:18;
y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2)+.05.*randn(size(x));
[Results,FitError]=peakfit([x;y],0,0,4,1,0,1,[4 2 9 2 12 2 14 2],0,0)
% 
% Sometimes experimental peaks are effected by exponential broadening,
% which does not by itself change the true peak areas, but does shift peak
% positions and increases peak width, overlap, and asymmetry, making it
% harder to separate the peaks. Fitting with a Gaussian does not work well:
%
y1=ExpBroaden(y',-15);
[Results,FitError]=peakfit([x;y1'],0,0,4,1,0,1,[4 2 9 2 12 2 14 2],0,0);
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on

Results
FitError
% Peakfit.m (and ipf.m) have an exponentially-broadened Gaussian peak shape
% (shape #5) that works well in those cases:, recovering the original peak
% positions, heights, and widths:
%
y1=ExpBroaden(y',-15);
disp('Calculating........')
[Results,FitError]=peakfit([x;y1'],0,0,4,5,15,1,[4 2 9 2 12 2 14 2],0,0);
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
Results
FitError
%---------------------------------------------------------------------
% Example 15: Prints out a table of parameter error estimates; Version 3
%  only. See DemoPeakfitBootstrap for a self-contained demo of this
%  function.
%
x=0:.05:9;
y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.01*randn(1,length(x));
[Results,FitError,Baseline,BestStart,xi,yi,BootstrapErrors]=peakfit([x;y],0,0,2,6,0,1,0,0,0);
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 16: (Version 3.2 or later)%  Fits both peaks of the Humps
%  function with a Gaussian/Lorentzian blend  (shape 13) 
%  that is 20% Gaussian. The 'Extra' argument sets the 
%  percentage of Gaussian shape (Extra=20).
%
x=0:.005:1;
y=humps(x);
[Results,FitError]=peakfit([x' y'],0,0,2,13,20,1,0,0)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 17:  (Version 3.2 or later)
%  Fit a slightly asymmetrical peak with a bifurcated Gaussian (shape 14).
%  The 'Extra' argument (=45) controls the peak asymmetry (50 is
%  symmetrical).
%
x=[0:.1:10];y=exp(-(x-4).^2)+.5*exp(-(x-5).^2)+.01*randn(size(x));
[Results,FitError]=peakfit([x' y'],0,0,1,14,45,10,0,0,0) 
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 18:  (Version 3.3 or later)
%  Example 1 without plotting or command window printing (11th input
%  argument = 0, default is 1)
% 
clf
x=[0:.1:10]';
y=exp(-(x-5).^2);peakfit([x y],0,0,1,1,0,0,0,0,0,0)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 19:  (Version 3.9 or later)
% Exponentially broadened Lorentzian with position=9, height=1,
% exponential factor ('extra') = 10.
x=0:.1:20; 
L=lorentzian(x,9,1);
L1=ExpBroaden(L',-10)+0.02.*randn(size(x))';
[Results,FitError]=peakfit([x;L1'],0,0,1,18,10,[9 1])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 20:  (Version 3.6 or later)
% Fixed-position Gaussian (shape 16), positions=[3 5]. 
x=0:.1:10;y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
[FitResults,FitError]=peakfit([x' y'],0,0,2,16,0,0,0,0,[3 5])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 21:  (Version 3.9 or later) Fitting the humps function with two
% Voigt profiles, flat baseline mode
[FitResults,FitError]=peakfit(humps(0:.01:2),71,140,2,20,1.7,1,[31 4.7 90 8.8],3)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 22: (Version 4.3 or later) Set +/- mode to 1 (bipolar)
x=0:.1:10;y=exp(-(x-5).^2)-.5*exp(-(x-3).^2)+.1*randn(size(x));
peakfit([x' y'],0,0,2,1,0,1,0,0,0,1,1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 23: Version 5 or later. Fits humps function to a model consisting 
% of one Pearson (shape=4, extra=3) and one Gaussian (shape=1), flat
% baseline mode=3, NumTrials=10.
x=0:.005:1.2;y=humps(x);
[FitResults,FitError]=peakfit([x' y'],0,0,2,[2 1],[0 0])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 24:  5 peaks, 5 different shapes, all heights = 1, widths = 3.
x=0:.1:60;
y=modelpeaks2(x,[1 2 3 4 5],[1 1 1 1 1],[10 20 30 40 50],[3 3 3 3 3],[0 0 0 2 -20])+.01*randn(size(x));
[FitResults,FitError]=peakfit([x' y'],0,0,5,[1 2 3 4 5],[0 0 0 2 -20])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 25:  Minimum width constraint (13th input argument)
x=1:30;y=gaussian(x,15,8)+.05*randn(size(x));
% No constraint:
[FitResults,FitError]=peakfit([x;y],0,0,5,1,0,10,0,0,0,1,0,0)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
% Widths constrained to values above 7:
[FitResults,FitError]=peakfit([x;y],0,0,5,1,0,10,0,0,0,1,0,7)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 26: Noise test with peak height = RMS noise = 1.
x=[-5:.02:5];y=exp(-(x).^2)+randn(size(x));P=peakfit([x;y],0,0,1,1,0,10,0,0,0,1,1);
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 27: Gaussian peak on strong sloped straight-line baseline, 2-peak
% fit with variable-slope straight line (shape 26, peakfit version 6 only).
x=8:.05:12;y=x+exp(-(x-10).^2);
[FitResults,FitError]=peakfit([x;y],0,0,2,[1 26],[1 1],1,0)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 28: Segmented linear fit (Shape 29, peakfit version 6 only):
x=[0.9:.005:1.7];y=humps(x);
peakfit([x' y'],0,0,9,29,0,10,0,0,0,1,1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 29: Polynomial fit (Shape 27, peakfit version 6 only)
x=[0.3:.005:1.7];y=humps(x);y=y+cumsum(y);
peakfit([x' y'],0,0,4,1,6,10,0,0,0,1,1)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 30: Effect of quantization of independent (x) and dependent (y) variables.
x=.5+round([2:.02:7.5]);y=exp(-(x-5).^2)+randn(size(x))/10;peakfit([x;y])
x=[2:.01:8];y=exp(-(x-5).^2)+.1.*randn(size(x));y=.1.*round(10.*y);peakfit([x;y])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 31: Variable-alpha Voigt functions, shape 30. (Version 7.45 and
% above only). FitResults has an added 6th column for the measured alphas
% of each peak.
x=[0:.005:1.3];y=humps(x);[FitResults,FitError]=peakfit([x' y'],60,120,2,30,15)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 32: Variable time constant exponentially braodened Gaussian
% functions, shape 31. (Version 7 and above only). FitResults has an added
% 6th column for the measured  time constant of each peak.
[FitResults,FitError]=peakfit(DataMatrix3,1860.5,364,2,31,32.9731,5,[1810 60 30 1910 60 30])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 33 Pearson variable shape, PeakShape=32,(Version 7 and above
% only). Requires modelpeaks2 function in path.
x=1:.1:30;y=modelpeaks2(x,[4 4],[1 1],[10 20],[5 5],[1 10]);
[FitResults,FitError]=peakfit([x;y],0,0,2,32,10,5)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 34 Gaussian/Lorentzian blend variable shape, PeakShape=33
% (Version 7 and above only). Requires modelpeaks2 function in path.
x=1:.1:30;y=modelpeaks2(x,[13 13],[1 1],[10 20],[3 3],[20 80]);
[FitResults,FitError]=peakfit([x;y],0,0,2,33,0,5)
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 35:  Fixed-position Gaussian (shape 16), positions=[3 5]. 
x=0:.1:10;y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
[FitResults,FitError]=peakfit([x' y'],0,0,2,16,0,0,0,0,[3 5])
echo off
fprintf(2,'Press any key to continue to next example.....or press Ctrl-C to abort.\n')
drawnow;pause
disp(' ')
echo on
%---------------------------------------------------------------------
% Example 36: Fixed-width Gaussian/Lorentzian blend (shape 35), 
% width=[3 3], compared to variable-width fit (shape 13). Fitting
% error is larger but nevertheless results are more accurate.
x=0:.1:10;y=GL(x,4,3,50)+.5*GL(x,6,3,50)+.1*randn(size(x));
[FitResults,FitError]=peakfit([x;y],0,0,2,13,50,1)
disp(' ')
[FitResults,FitError]=peakfit([x;y],0,0,2,35,50,1,0,0,[3 3])
echo off
drawnow
disp(' ')
echo 
%---------------------------------------------------------------------
echo off
disp('*********** End of testpeakfit.m ****************')
