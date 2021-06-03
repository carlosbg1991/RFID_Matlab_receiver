function PowerTransformCalibrationCurve(power)
% PowerMethodCalibrationCurve simulates the calibration of a flow injection
% system or other technique that produces signal peaks that are related to
% an underlying concentration or amplitude ('amp'). In this example, six
% known standards are measured sequentially, resulting in six peaks in the
% observed signal. To simulate a realistic measurement, three sources of
% random disturbance are added to the observed signal: Noise - random white
% noise added to the signal data points; background - broad curved random
% background; and broadening - random exponential broadening of the peaks.
% In the power method, the raw time-series signal is raised to a positive
% power ('power' in line 42) in order to sharpen the peaks and reduce overlap. The
% resulting peaks are measured by the measurepeaks.m function, which
% determines the absolute peak height, the peak-valley difference, the
% perperdicular drop area, and the tangent skim area. Each of these
% measures are raised to the 1/power in order to compensate for the power,
% then compared by plotting them in figures 2-5 against the true underlying
% amplitudes ("amp") and computing the slope, intercept, and R2. If all the
% random disturbances are set to zero (in lines 33-36), the R2 values will
% all be 1.000. Otherwise the measurements will not be perfect and some
% methods will result in better measurements (R2 closer to 1.000) than
% others. Note: you can expand the number of calibration standards by
% adding additinal elements to "amp", "pos", and "wid". 
% You must have gaussian.m, expgaussian.m, measurepeaks.m, and plotit.m in
% the path. Download from http://tinyurl.com/cey8rwh
% Tom O'Haver (toh@umd.edu) Version 2, September 2018
%
format compact
format SHORT G
increment=1; % Interval between data points in the signal
x=1:increment:700; % defines the "time" axis.

% For each simulated peak, compute the amplitude, position, and width
% amp=[1 2 3 4 5 6];  % True amplitudes of the peaks  (Change if desired)
amp=2+round(6.*rand(1,6));  % True amplitudes of the peaks are random
pos=[100 200 300 400 500 600];   % Positions of the peaks (Change if desired)
wid=70.*ones(size(amp));   % Widths of the peaks (Change if desired)
% Signal disturbances. Change as desired: 0=perfect situation; larger=bad
Noise=0.0; % Amount of random noise added to the signal.  
background=0; % Amplitude of broad random background/baseline under the peaks
broadening=0; % Random broadening of peaks; try 0 to 5  
FinalSmooth=6; % Final smoothing by the detector or electronics 
% power=2; % Exponent for power method (1-normal method)
% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      if broadening>0,
          A(k,:)=expgaussian(x,pos(k),wid(k),-(2+broadening*rand(1)))+.5.*expgaussian(x,pos(k)-10,wid(k),-(2+broadening*rand(1)));
      else
          A(k,:)=gaussian(x,pos(k),wid(k));
      end
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)];
      TA(p,:)=trapz(x,A(k,:));
      p=p+1;
end
TrueAreas=amp.*TA';
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
nz=z+Noise.*randn(size(z));  % Optionally adds random noise
ry=nz+background*gaussian(x,700*rand(1),800); % Adds a broad background whose amplitude=background
if FinalSmooth<1;FinalSmooth=1;end
ys=fastsmooth(ry,FinalSmooth,3,1); % Final smoothing
y=(ys.^power);  % Power method
y=y.*max(ys)./max(y);

% Detect and measure peaks in the signal
SlopeThreshold=.000001; % 
AmpThreshold=.01;
SmoothWidth=7;  % SmoothWidth should be roughly equal the peak width (in points)
FitWidth=11; % FitWidth should be roughly equal to the peak widths (in points)
sizex=size(x);
NumPoints=sizex(2);
figure(1)
Measured=measurepeaks(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,0);
sizeP=size(Measured);
NumPeaks=sizeP(1);
plot(x,ys,x,y)
disp(' ')
disp(['PowerMethodCalibrationCurve.m    Power = ' num2str(power) ])
title('Power transofrm method demo. Blue = Original signal.   Green = Power transform')
% Label the x-axis with the parameter values
xlabel(['Power=' num2str(power) '    SlopeThresh.=' num2str(SlopeThreshold) '    AmpThresh.=' num2str(AmpThreshold) '    SmoothWidth=' num2str(SmoothWidth) '    FitWidth=' num2str(FitWidth)])
clear PeakHeight
clear PeakTable

disp('          Peak    Position      PeakMax    Peak-valley    Perp drop    Tan skim')

disp(Measured)
figure(2)
clf
[coef,APM]=plotit(amp,Measured(:,3).^(1/power),1);
title(['Absolute Peak Maximum.  Power=' num2str(power) ])
ylabel('Measured amplitude')
xlabel('True amplitude')

figure(3)
clf
[coef,PVD]=plotit(amp,Measured(:,4).^(1/power),1);
title(['Peak - Valley difference      Power=' num2str(power) ])
xlabel('True amplitude')
ylabel('Measured amplitude')

figure(4)
clf
[coef,PDA]=plotit(amp,Measured(:,5).^(1/power),1);
title(['Perpendicular drop Area.     Power=' num2str(power) ])
xlabel('True amplitude')
ylabel('Measured peak area')

figure(5)
clf
[coef,TSA]=plotit(amp,Measured(:,6).^(1/power),1);
title(['Tangent skim Area.     Power=' num2str(power) ])
xlabel('True amplitude')
ylabel('Measured peak area')

figure(6)
clf
plotit(TrueAreas,Measured(:,6).^(1/power),1);
title(['Correlation between true and measured tangent skim area.    Power=' num2str(power) ])
xlabel('True peak area')
ylabel('Measured peak area')

disp('R2 values')
disp('      PeakMax    Peak-valley    Perp drop    Tan skim')
disp([APM PVD PDA TSA])

function [FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight)
% A command-line peak fitting program for time-series signals, written as a
% self-contained Matlab function in a single m-file. Uses a non-linear
% optimization algorithm to decompose a complex, overlapping-peak signal
% into its component parts. The objective is to determine whether your
% signal can be represented as the sum of fundamental underlying peaks
% shapes. Accepts signals of any length, including those with non-integer
% and non-uniform x-values. Fits any number of peaks of any of 45 curve
% shapes. This is a command line version, usable from a remote terminal. It
% is capable of making multiple trial fits with sightly different starting
% values and taking the one with the lowest mean fit error (example 6). It
% can estimate the standard deviation of peak parameters from a single
% signal using the bootstrap method (example 10).
%
% If you are unsure what input arguments to use (number of peaks, shape,
% etc.), try fitting your signal using the interactive peak fitter ipf.m
% (Matlab only), which uses single keystrokes to pan and zoom, select
% number of peaks, peak shape, baseline mode, etc. Then press W to print
% out the peakfit command line for that fit.
%
% Important: the data matrix "signal" must be a 2xn or nx2 matrix; the x
% and y vectors must be separate rows or columns and not concatenated into
% one long vector.
%
% Version 8.3: February 2017, correction to BiGaussian shape.
%
% For more details, see
% http://terpconnect.umd.edu/~toh/spectrum/CurveFittingC.html and
% http://terpconnect.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
%
% peakfit(signal);       
% Performs an iterative least-squares fit of a single Gaussian peak to the
% data matrix "signal", which has x values in column 1 and Y values in
% column 2 (e.g. [x y]). 
%
% peakfit(signal,center,window);
% Fits a single Gaussian peak to a portion of the matrix "signal". The
% portion is centered on the x-value "center" and has width "window" (in x
% units).
% 
% peakfit(signal,center,window,NumPeaks..);
% "NumPeaks" = number of peaks in the model (default is 1 if not
% specified). No limit to maximum number of peaks in version 3.1
% 
% peakfit(signal,center,window,NumPeaks,peakshape); "peakshape" specifies
% the peak shape of the model: (1=Gaussian (default), 2=Lorentzian,
% 3=logistic distribution, 4=Pearson, 5=exponentionally broadened Gaussian;
% 6=equal-width Gaussians; 7=Equal-width Lorentzians; 8=exponentionally
% broadened equal-width Gaussian, 9=exponential pulse, 10=up-sigmoid
% (logistic function), 11=Fixed-width Gaussian, 12=Fixed- width Lorentzian;
% 13=Gaussian/Lorentzian blend; 14=Bifurcated Gaussian,
% 15=Breit-Wigner-Fano, 16=Fixed-position Gaussians; 17=Fixed-position 
% Lorentzians; 18=exponentionally broadened Lorentzian; 19=alpha function; 
% 20=Voigt profile; 21=triangular; 23=down-sigmoid; 25=lognormal; 26=slope;
% 27=Gaussian first derivative; 28=polynomial; 29=piecewise linear;
% 30=variable-alpha Voigt; 31=variable time constant ExpGaussian;
% 32=variable shape Pearson; 33=variable Gaussian/Lorentzian blend,
% 34=fixed-width Voigt; 35=fixed-width Gaussian/Lorentzian blend; 
% 36=fixed-width exponentionally-broadened Gaussian; 37=fixed-width 
% Pearson; 38=variable time constant ExpLorentzian; 40=sine wave;
% 41=rectangle; 42=flattened Gaussian; 43:Gompertz function (3 variables);
% 44=1-exp(-k*x); 45: Four-parameter logistic; 46=quadratic/parabola;
% 47=blackbody emission; 48=equal-width exponential pulses.
%
% peakfit(signal,center,window,NumPeaks,peakshape,extra) 
% 'extra' specifies the value of 'extra', used only in the Voigt, Pearson,
% exponentionally broadened Gaussian, Gaussian/Lorentzian blend, and
% bifurcated Gaussian and Lorentzian shapes to fine-tune the peak shape.
% 
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials);
% Performs "NumTrials" trial fits and selects the best one (with lowest
% fitting error). NumTrials can be any positive integer (default is 1).
%
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start)
% Specifies the first guesses vector "firstguess" for the peak positions
% and widths. Must be expressed as a vector , in square brackets, e.g.
% start=[position1 width1 position2 width2 ...]
% 
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero) 
% 'autozero' sets the baseline correction mode:
% autozero=0 (default) does not subtract baseline from data segment; 
% autozero=1 interpolates a linear baseline from the edges of the data 
%            segment and subtracts it from the signal (assumes that the 
%            peak returns to the baseline at the edges of the signal); 
% autozero=2 is like mode 1 except that it computes a quadratic curved baseline; 
% autozero=3 compensates for a flat baseline without reference to the 
%            signal itself (best if the peak does not return to the
%            baseline at the edges of the signal).
%
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
% autozero,fixedparameters)
% 'fixedparameters' specifies fixed values for widths (shapes 10, 11) or
% positions (shapes 16, 17)
% 
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
% autozero,fixedparameters,plots)
% 'plots' controls graphic plotting: 0=no plot; 1=plots draw as usual (default)
% 
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
% autozero,fixedparameters,plots,bipolar)
% 'bipolar' = 0 constrain peaks heights to be positions; 'bipolar' = 1
% allows positive ands negative peak heights.
%
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
% autozero,fixedparameters,plots,bipolar,minwidth)
% 'minwidth' sets the minmimum allowed peak width. The default if not
% specified is equal to the x-axis interval. Must be a vector of minimum
% widths, one value for each peak, if the multiple peak shape is chosen, 
% as in example 17 and 18.
%
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
% autozero,fixedparameters,plots,bipolar,minwidth)
%  'DELTA' (14th input argument) controls the restart variance when
%  NumTrials>1. Default value is 1.0. Larger values give more variance.
%  Version 5.8 and later only. 
% [FitResults,FitError]=peakfit(signal,center,window...) Returns the
% FitResults vector in the order peak number, peak position, peak height,
% peak width, and peak area), and the FitError (the percent RMS
% difference between the data and the model in the selected segment of that
% data) of the best fit.
%
% [FitResults,LowestError,BestStart,xi,yi,BootResults]=peakfit(signal,...)
% Prints out parameter error estimates for each peak (bootstrap method).
%
% Optional output parameters 
% 1. FitResults: a table of model peak parameters, one row for each peak,
%    listing Peak number, Peak position, Height, Width, and Peak area.
% 2. GOF: Goodness of Fit, a vector containing the rms fitting error of the
%    best trial fit and the R-squared (coefficient of determination).
% 3. Baseline, the polynomial coefficients of  the baseline in linear
%    and quadratic baseline modes (1 and 2) or the value of the constant
%    baseline in flat baseline mode.
% 3. coeff: Coefficients for the polynomial fit (shape 28 only; for other 
%    shapes, coeff=0)
% 5. residual: the difference between the data and the best fit.
% 6. xi: vector containing 600 interploated x-values for the model peaks. 
% 7. yi: matrix containing the y values of each model peak at each xi. 
%    Type plot(xi,yi(1,:)) to plot peak 1 or plot(xi,yi) to plot all peaks
% 8. BootResults: a table of bootstrap precision results for a each peak
%    and peak parameter.
%  
% Example 1a: Fits exp(-x)^2 with a single Gaussian peak model.
% >> x=[0:.1:10]';y=exp(-(x-5).^2);peakfit([x y])
%    Peak number  Peak position   Height     Width      Peak area
%         1            5            1        1.665       1.7725
%
% Example 1b: Fits histogram of 10000 random numbers with a single Gaussian
% >> [N,X]=hist(randn(size(1:10000)));peakfit([X;N])
%   Peak number position  Height     Width     Area
%     1.0000    0.0659   234.4104    2.4592   607.3360
%
% Example 1c: Fits small set of manually entered y data to a single Gaussian peak model.
% >> y=[0 1 2 4 6 7 6 4 2 1 0 ];x=1:length(y);
% >> peakfit([x;y],length(y)/2,length(y),0,0,0,0,0,0)
%
% Example 2:
% x=[0:.01:10];y=exp(-(x-5).^2)+randn(size(x));peakfit([x;y])
% Measurement of very noisy peak with signal-to-noise ratio = 1.
% ans =
%         1         5.0279       0.9272      1.7948      1.7716
%
% Example 3:
% >> x=[0:.1:10];y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
% >> peakfit([x' y'],0,0,2)
% Fits a noisy two-peak signal with a double Gaussian model (NumPeaks=2).
% ans =
%         1       3.0001      0.49489        1.642      0.86504
%         2       4.9927       1.0016       1.6597       1.7696
%
% Example 4:
% >> x=1:10;y=ones(size(x))./(1+(x-5).^2);peakfit(y,0,0,1,2)
% Fit Lorentzian (peakshape=2) located at x=50, height=1, width=2.
% ans =
%        1           50      0.99974       1.9971       3.1079
%
% Example 5: 
% >> x=[0:.005:1];y=humps(x);peakfit([x' y'],.3,.7,1,4,3);
% Fits a portion of the humps function, 0.7 units wide and centered on 
% x=0.3, with a single (NumPeaks=1) Pearson function (peakshape=4) with
% extra=3 (controls shape of Pearson function).
%
% Example 6: 
% >> x=[0:.005:1];y=(humps(x)+humps(x-.13)).^3;smatrix=[x' y'];
% >> [FitResults,FitError]=peakfit(smatrix,.4,.7,2,1,0,10)
% Creates a data matrix 'smatrix', fits a portion to a two-peak Gaussian
% model, takes the best of 10 trials.  Returns FitResults and FitError.
% FitResults =
%          1      0.31056  2.0125e+006      0.11057  2.3689e+005
%          2      0.41529  2.2403e+006      0.12033  2.8696e+005
% FitError =
%         1.1899
%
% Example 7:
% >> peakfit([x' y'],.4,.7,2,1,0,10,[.3 .1 .5 .1]);
% As above, but specifies the first-guess position and width of the two
% peaks, in the order [position1 width1 position2 width2]
%
% Example 8: (Version 4 only)
% Demonstration of the four autozero modes, for a single Gaussian on flat
%  baseline, with position=10, height=1, and width=1.66. Autozero mode
%  is specified by the 9th input argument (0,1,2, or 3).
% >> x=8:.05:12;y=1+exp(-(x-10).^2);
% >> [FitResults,FitError]=peakfit([x;y],0,0,1,1,0,1,0,0)
% Autozero=0 means to ignore the baseline (default mode if not specified)
% FitResults =
%         1           10       1.8561        3.612       5.7641
% FitError =
%         5.387
% >> [FitResults,FitError]=peakfit([x;y],0,0,1,1,0,1,0,1)
% Autozero=1 subtracts linear baseline from edge to edge.
% Does not work well because signal does not return to baseline at edges.
% FitResults =
%         1       9.9984      0.96153        1.559       1.5916
% FitError =
%        1.9801
% >> [FitResults,FitError]=peakfit([x;y],0,0,1,1,0,1,0,2)
% Autozero=1 subtracts quadratic baseline from edge to edge.
% Does not work well because signal does not return to baseline at edges.
% FitResults =
%         1       9.9996      0.81749       1.4384       1.2503
% FitError =
%        1.8204
% Autozero=3: Flat baseline mode, measures baseline by regression
% >> [FitResults,FitError,Baseline]=peakfit([x;y],0,0,1,1,0,1,0,3)
% FitResults =
%         1           10       1.0001       1.6653       1.7645
% Baseline =
%     0.0037056
% FitError =
%       0.99985
%
% Example 9:
% x=[0:.1:10];y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
% [FitResults,FitError]=peakfit([x' y'],0,0,2,11,0,0,0,0,[1.666 1.666])
% Same as example 3, fit with fixed-width Gaussian (shape 11), width=1.666
% 
% Example 10: (Version 3 or later; Prints out parameter error estimates)
% x=0:.05:9;y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.01*randn(1,length(x));
% [FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit([x;y],0,0,2,6,0,1,0,0,0);
%
% Example 11: (Version 3.2 or later)
% x=[0:.005:1];y=humps(x);[FitResults,FitError]=peakfit([x' y'],0.54,0.93,2,13,15,10,0,0,0) 
%
% FitResults =
%         1      0.30078       190.41      0.19131       23.064
%         2      0.89788       39.552      0.33448       6.1999
% FitError = 
%       0.34502
% Fits both peaks of the Humps function with a Gaussian/Lorentzian blend
% (shape 13) that is 15% Gaussian (Extra=15).
% 
% Example 12:  (Version 3.2 or later)
% >> x=[0:.1:10];y=exp(-(x-4).^2)+.5*exp(-(x-5).^2)+.01*randn(size(x));
% >> [FitResults,FitError]=peakfit([x' y'],0,0,1,14,45,10,0,0,0) 
% FitResults =
%         1       4.2028       1.2315        4.077       2.6723
% FitError =
%       0.84461
% Fit a slightly asymmetrical peak with a bifurcated Gaussian (shape 14)
% 
% Example 13:  (Version 3.3 or later)
% >> x=[0:.1:10]';y=exp(-(x-5).^2);peakfit([x y],0,0,1,1,0,0,0,0,0,0)
% Example 1 without plotting (11th input argument = 0, default is 1)
% 
% Example 14:  (Version 3.9 or later)
% Exponentially broadened Lorentzian with position=9, height=1.
% x=[0:.1:20]; 
% L=lorentzian(x,9,1);
% L1=ExpBroaden(L',-10)+0.02.*randn(size(x))';
% [FitResults,FitError]=peakfit([x;L1'],0,0,1,18,10)
%
% Example 15: Fitting humps function with two Voigts (version 7.45)
% FitResults,FitError]=peakfit(humps(0:.01:2),60,120,2,20,1.7,5,0)
% FitResults =
%         1       31.099       94.768       19.432       2375.3
%         2       90.128       20.052       32.973       783.34
% FitError =
%       0.72829      0.99915
%
% Example 16: (Version 4.3 or later) Set +/- mode to 1 (bipolar)
% x=[0:.1:10];y=exp(-(x-5).^2)-.5*exp(-(x-3).^2)+.1*randn(size(x));
% peakfit([x' y'],0,0,2,1,0,1,0,0,0,1,1)
% FitResults =
%         1       3.1636      -0.5433         1.62      -0.9369
%         2       4.9487      0.96859       1.8456       1.9029
% FitError =
%        8.2757
%
% Example 17: Version 5 or later. Fits humps function to a model consisting 
% of one Lorentzian (shape=2) and one Gaussian (shape=1). 
% x=[0:.005:1.2];y=humps(x);[FitResults,FitError]=peakfit([x' y'],0,0,2,[2 1],[0 0])
% FitResults =
%     1.0000    0.3018   97.3992    0.1886   25.1164
%     2.0000    0.8962   18.7892    0.3369    6.6241
% FitError =
%     1.0744    0.998
%
% Example 18: 5 peaks, 5 different shapes, all heights = 1, widths = 3.
% x=0:.1:60;
% y=modelpeaks2(x,[1 2 3 4 5],[1 1 1 1 1],[10 20 30 40 50],...
% [3 3 3 3 3],[0 0 0 2 -20])+.01*randn(size(x));
% peakfit([x' y'],0,0,5,[1 2 3 4 5],[0 0 0 2 -20])
%
% Example 19: Minimum width constraint (13th input argument)
% x=1:30;y=gaussian(x,15,8)+.05*randn(size(x));
% No constraint:
% peakfit([x;y],0,0,5,1,0,10,0,0,0,1,0,0);
% Widths constrained to values above 7:
% peakfit([x;y],0,0,5,1,0,10,0,0,0,1,0,7);
%
% Example 20: Noise test with peak height = RMS noise = 1.
% x=[-5:.02:5];y=exp(-(x).^2)+randn(size(x));P=peakfit([x;y],0,0,1,1,0,10,0,0,0,1,1);
%
% Example 21a: Gaussian peak on strong sloped straight-line baseline, 2-peak
% fit with variable-slope straight line ("third peak" is shape 26, peakfit version 6 only).
% x=8:.05:12;y=x+exp(-(x-10).^2);
% [FitResults,FitError]=peakfit([x;y],0,0,2,[1 26],[1 1],1,0)
% FitResults =           
%         1           10            1       1.6651       1.7642
%         2        4.485      0.22297         0.05       40.045
% FitError =0.093
%
% Example 21b: As above, but with TWO signal peaks, peakshape=[1 1 26],
% helped by adding rough first guesses {'start') using the 'polyfit'
% function to generate first guesses for the sloping baseline:
%  x=7:.05:13;y=x/2+exp(-(x-9).^2)+exp(-(x-11).^2);
% start=[8 1 10 1 polyfit(x,y,1)];
% [FitResults,FitError]=peakfit([x;y],0,0,3,[1 1 26],[1 1 1],1,start)
%
% Example 22: Segmented linear fit (Shape 29, peakfit version 6 only):
% x=[0.9:.005:1.7];y=humps(x);
% peakfit([x' y'],0,0,9,29,0,10,0,0,0,1,1)
%
% Example 23: Sixth degree polynomial fit (Shape 28, peakfit version 6 only)
%   x=[0.3:.005:1.7];y=humps(x);y=y+cumsum(y);
%   peakfit([x' y'],0,0,1,28,6,10,0,0,0,1,1)
%
% Example 24: Effect of quantization of independent (x) and dependent (y) variables.
% x=.5+round([2:.02:7.5]);y=exp(-(x-5).^2)+randn(size(x))/10;peakfit([x;y])
% x=[2:.01:8];y=exp(-(x-5).^2)+.1.*randn(size(x));y=.1.*round(10.*y);peakfit([x;y])
%
% Example 25: Variable-alpha Voigt functions, shape 30. (Version 7.45 and
% above only). FitResults has an added 6th column for the measured alphas
% of each peak.
% x=[0:.005:1.3];y=humps(x);[FitResults,FitError]=peakfit([x' y'],60,120,2,30,15)
% FitResults =
%      1      0.30147       95.797      0.19391       24.486       2.2834
%      2      0.89186       19.474      0.33404       7.3552      0.39824
% FitError =
%       0.84006      0.99887
%
% Example 26: Variable time constant exponentially-broadened Gaussian
% function, shape 31. (Version 7 and above only). FitResults has an added
% 6th column for the measured  time constant of each peak.
% [FitResults,FitError]=peakfit(DataMatrix3,1860.5,364,2,31,32.9731,5,[1810 60 30 1910 60 30])
% FitResults =
%     1       1800.6       1.8581       60.169       119.01       32.781
%     2       32.781      0.48491       1900.4        30.79       33.443
% FitError =
%      0.076651      0.99999
%
% Example 27 Pearson variable shape, PeakShape=32,(Version 7 and above
% only). Requires modelpeaks2 function in path.
% x=1:.1:30;y=modelpeaks2(x,[4 4],[1 1],[10 20],[5 5],[1 10]);
% [FitResults,FitError]=peakfit([x;y],0,0,2,32,10,5)
%
% Example 28 Gaussian/Lorentzian blend variable shape, PeakShape=33
% (Version 7 and above only). Requires modelpeaks2 function in path.
% x=1:.1:30;y=modelpeaks2(x,[13 13],[1 1],[10 20],[3 3],[20 80]);
% [FitResults,FitError]=peakfit([x;y],0,0,2,33,0,5)
%
% Example 29:  Fixed-position Gaussian (shape 16), positions=[3 5]. 
% x=0:.1:10;y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
% [FitResults,FitError]=peakfit([x' y'],0,0,2,16,0,0,0,0,[3 5])
%
% Example 30: Fixed-width Gaussian/Lorentzian blend (shape 35), 
% width=[3 3], compared to variable-width fit (shape 13). Fitting
% error is larger but nevertheless results are more accurate.
% x=0:.1:10;y=GL(x,4,3,50)+.5*GL(x,6,3,50)+.1*randn(size(x));
% [FitResults,FitError]=peakfit([x;y],0,0,2,13,50,1)
%          1       4.0632       1.0545       3.2182       3.9242
%          2       6.2736      0.41234       2.8114       1.3585
%  FitError =  6.4311      0.95211 
% [FitResults,FitError]=peakfit([x;y],0,0,2,35,50,1,0,0,[3 3])
% FitResults =
%          1       3.9527       1.0048            3       3.5138
%          2       6.1007       0.5008            3       1.7502
% 
%  FitError =  6.4783      0.95141
%
% Example 31: Variable time constant exponentially-broadened Lorentzian
% function, shape 38. (Version 7.7 and above only). FitResults has an added
% 6th column for the measured time constant.
% x=[1:100]';
% y=explorentzian(x',40,5,-10)+.01*randn(size(x));
% peakfit([x y],0,0,1,38,0,10)
%
% Example 32: 3-parameter logistic (Gompertz), shape 43. (Version 7.9 and
% above only). Parameters labeled Bo, Kh, and L. FitResults extended to 6
% columns. 
% t=0:.1:10; Bo=6;Kh=3;L=4;
% y=Bo*exp(-exp((Kh*exp(1)/Bo)*(L-t)+1))+.1.*randn(size(t));
% [FitResults,GOF]=peakfit([t;y],0,0,1,43)
%
% Example 33: Shape 46, 'quadslope' (version 8 and above). Two overlapping
% Gaussians (positions=9,11; heights=1; widths=1.666) on a curved baseline,
% using a 3-peak fit with peakshape=[1 1 46], default NumTrials and start.
% x=7:.05:13;
% y=x.^2/50+exp(-(x-9).^2)+exp(-(x-11).^2)+.01.*randn(size(x));
% [FitResults,FitError]=peakfit([x;y],0,0,3,[1 1 46],[1 1 1])

% Copyright (c) 2015, Thomas C. O'Haver
 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
% 
global AA xxx PEAKHEIGHTS FIXEDPARAMETERS AUTOZERO delta BIPOLAR CLIPHEIGHT
% format short g
format compact
warning off all
NumArgOut=nargout;
datasize=size(signal);
if datasize(1)<datasize(2),signal=signal';end
datasize=size(signal);
if datasize(2)==1, %  'signal' is a vector; Must be peakfit(Y-vector)
    X=1:length(signal); % Create an independent variable vector
    Y=signal;
else
    % 'signal' is a matrix. Must be peakfit(DataMatrix)
    X=signal(:,1); % Split matrix argument 
    Y=signal(:,2);
end
X=reshape(X,1,length(X)); % Adjust X and Y vector shape to 1 x n (rather than n x 1)
Y=reshape(Y,1,length(Y));
% If necessary, flip the data vectors so that X increases
if X(1)>X(length(X)),
    disp('X-axis flipped.')
    X=fliplr(X);
    Y=fliplr(Y);
end

% Isolate desired segment from data set for curve fitting
if nargin==1 || nargin==2,center=(max(X)-min(X))/2;window=max(X)-min(X);end
% Y=Y-min(Y);
xoffset=0;
n1=val2ind(X,center-window/2);
n2=val2ind(X,center+window/2);
if window==0,n1=1;n2=length(X);end
xx=X(n1:n2)-xoffset;
yy=Y(n1:n2);
ShapeString='Gaussian';
coeff=0;
CLIPHEIGHT=max(Y);
LOGPLOT=0;
% Define values of any missing arguments
% (signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA)
switch nargin
    case 1
        NumPeaks=1;
        peakshape=1;
        extra=1;
        NumTrials=1;
        xx=X;yy=Y;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 2
        NumPeaks=1;
        peakshape=1;
        extra=1;
        NumTrials=1;
        xx=signal;yy=center;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 3
        NumPeaks=1;
        peakshape=1;
        extra=1;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 4 % Numpeaks specified in arguments
        peakshape=1;
        extra=1;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 5 % Numpeaks, peakshape specified in arguments
        extra=ones(1,NumPeaks);
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 6
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 7
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 8
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 9
        % initialstart=start % testing
        AUTOZERO=autozero;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 10
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 11
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 12
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 13
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=1;
    case 14
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=DELTA;
        CLIPHEIGHT=max(Y);
    case 15
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=DELTA;
        CLIPHEIGHT=clipheight;
    otherwise
end % switch nargin
% Saturation Code, skips points greater than set maximum
if CLIPHEIGHT<max(Y),
    apnt=1;
    for pnt=1:length(xx),
        if yy(pnt)<CLIPHEIGHT,
            axx(apnt)=xx(pnt);
            ayy(apnt)=yy(pnt);
            apnt=apnt+1;
        end
    end
    xx=axx;yy=ayy;
end
% Default values for placeholder zeros1
if NumTrials==0;NumTrials=1;end
shapesvector=peakshape;
if isscalar(peakshape),
else
    % disp('peakshape is vector');
    shapesvector=peakshape;
    NumPeaks=length(peakshape);
    peakshape=22;
end
if peakshape==0;peakshape=1;end
if NumPeaks==0;NumPeaks=1;end
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
firststart=start; % <<<<<<<<<<<
if FIXEDPARAMETERS==0, FIXEDPARAMETERS=length(xx)/10;end
if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
if peakshape==17;FIXEDPOSITIONS=fixedparameters;end
if AUTOZERO>3,AUTOZERO=3;disp('AUTPOZERO must be between 0 and 3');end
if AUTOZERO<0,AUTOZERO=0;disp('AUTPOZERO must be between 0 and 3');end
Heights=zeros(1,NumPeaks);
FitResults=zeros(NumPeaks,6);

% % Remove linear baseline from data segment if AUTOZERO==1
baseline=0;
bkgcoef=0;
bkgsize=round(length(xx)/10);
if bkgsize<2,bkgsize=2;end
lxx=length(xx);
if AUTOZERO==1, % linear autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:(round(length(xx)/bkgsize)));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1);  % Fit straight line to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if
if AUTOZERO==2, % Quadratic autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero

PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;
% Assign ShapStrings
switch peakshape(1)
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';
    case 3
        ShapeString='Logistic';
    case 4
        ShapeString='Pearson';
    case 5
        ShapeString='ExpGaussian';
    case 6
        ShapeString='Equal width Gaussians';
    case 7
        ShapeString='Equal width Lorentzians';
    case 8
        ShapeString='Exp. equal width Gaussians';
    case 9
        ShapeString='Exponential Pulse';
    case 10
        ShapeString='Up Sigmoid (logistic function)';
    case 23
        ShapeString='Down Sigmoid (logistic function)';  
    case 11
        ShapeString='Fixed-width Gaussian';
    case 12
        ShapeString='Fixed-width Lorentzian';
    case 13
        ShapeString='Gaussian/Lorentzian blend';
    case 14
        ShapeString='BiGaussian';    
    case 15
        ShapeString='Breit-Wigner-Fano';   
    case 16
        ShapeString='Fixed-position Gaussians';
    case 17
        ShapeString='Fixed-position Lorentzians';
    case 18
        ShapeString='Exp. Lorentzian';
    case 19
        ShapeString='Alpha function';
    case 20
        ShapeString='Voigt (equal alphas)';
    case 21
        ShapeString='triangular';
    case 22
        ShapeString=num2str(shapesvector);
    case 24
        ShapeString='Negative Binomial Distribution';
    case 25
        ShapeString='Lognormal Distribution';
    case 26
        ShapeString='slope';
    case 27
        ShapeString='First derivative';
    case 28
        ShapeString='Polynomial';
    case 29
        ShapeString='Segmented linear';
    case 30
        ShapeString='Voigt (variable alphas)';
    case 31
        ShapeString='ExpGaussian (var. time constant)';
    case 32
        ShapeString='Pearson (var. shape constant)';
    case 33
        ShapeString='Variable Gaussian/Lorentzian';
    case 34
        ShapeString='Fixed-width Voigt';
    case 35
        ShapeString='Fixed-width G/L blend';
    case 36
        ShapeString='Fixed-width ExpGaussian';
    case 37
        ShapeString='Fixed-width Pearson';
    case 38
        ShapeString='ExpLorentzian (var. time constant)'; 
    case 40
        ShapeString='Sine wave';
    case 41
        ShapeString='Rectangular pulse';
    case 42
        ShapeString='Flattened Gaussian';  
    case 43
        ShapeString='3-parameter Gompertz.';  
    case 44
        ShapeString='1-exp(-k*t)';  
    case 45
        ShapeString='4-parameter logistic'; 
    case 46
        ShapeString='Quadratic Baseline'; 
    case 47
        ShapeString='Blackbody'; 
    case 48
        ShapeString='Equal-width exp. pulse'; 
    otherwise
end % switch peakshape
  
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.0000001,'Display','off','MaxFunEvals',1000 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*3); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
TrialParameters=zeros(1,NumPeaks.*3);
for k=1:NumTrials,
    % StartMatrix(k,:)=newstart;
    % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
    switch peakshape(1)
        case 1
            TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 2
            TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 3
            TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 4
            TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 5
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 6
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 7
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewlorentzian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 8
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 9
            TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 10
            TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 23
            TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 11
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            % fixedstart=fixedstart
            TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart,options);
        case 12
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart,options);
        case 13
            TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 14
            TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 15
            TrialParameters=fminsearch(@(lambda)(fitBWF(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 16
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 17
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 18
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 19
            TrialParameters=fminsearch(@(lambda)(fitalphafunction(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 20
            TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 21
            TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 22
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,NumPeaks,shapesvector,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH(Peak),
                    TrialParameters(2*Peak)=MINWIDTH(Peak);
                end
            end
        case 24
            TrialParameters=fminsearch(@(lambda)(fitnbinpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 25
            TrialParameters=fminsearch(@(lambda)(fitlognpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 26
            TrialParameters=fminsearch(@(lambda)(fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options);
             coeff=TrialParameters;
        case 27
            TrialParameters=fminsearch(@(lambda)(fitd1gauss(lambda,xx,yy)),newstart,options);
        case 28
            coeff=fitpolynomial(xx,yy,extra);
            TrialParameters=coeff;
        case 29
            cnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cnewstart(pc)=newstart(2.*pc-1)+(delta*(rand-.5)/50);
            end
            TrialParameters=fminsearch(@(lambda)(fitsegmented(lambda,xx,yy)),cnewstart,options);
        case 30
            % newstart=newstart % testing
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
%             case30newstart=newstart % uncomment for testing
%             sizestartcase30=size(start) % uncomment for testing

            TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),start);
         case 31
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % case31newstart=newstart % uncomment for testing
            % sizestartcase31=size(start) % uncomment for testing
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussianv(lambda,zxx,zyy)),newstart);
        case 32
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            case32newstart=newstart % uncomment for testing
            sizestartcase32=size(start) % uncomment for testing
            TrialParameters=fminsearch(@(lambda)(fitpearsonv(lambda,xx,yy)),newstart);
        case 33
             nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart % uncomment for testing
            TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart);
        case 34
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end  
            % fixedstart=fixedstart % uncomment for testing
            TrialParameters=fminsearch(@(lambda)(fitFWVoigt(lambda,xx,yy,extra)),fixedstart,options);
        case 35
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWGL(lambda,xx,yy,extra)),fixedstart,options);
        case 36
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWExpGaussian(lambda,xx,yy,extra)),fixedstart,options);
        case 37
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWPearson(lambda,xx,yy,extra)),fixedstart,options);
        case 38
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzianv(lambda,zxx,zyy)),newstart);
        case 40
             TrialParameters=fminsearch(@(lambda)(fitsine(lambda,xx,yy)),newstart,options); 
        case 41
            TrialParameters=fminsearch(@(lambda)(fitrectangle(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 42
            TrialParameters=fminsearch(@(lambda)(fitngaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 43
             nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            TrialParameters=fminsearch(@(lambda)(fitGompertz(lambda,xx,yy)),newstart);

        case 44
            TrialParameters=fminsearch(@(lambda)(fitOneMinusExp(lambda,xx,yy)),newstart,options);
         case 45
             nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart % uncomment for testing
            TrialParameters=fminsearch(@(lambda)(fitFourPL(lambda,xx,yy)),newstart);
        case 46
            TrialParameters=fminsearch(@(lambda)(fitquadslope(lambda,xx,yy)),newstart,options);
        case 47
            bbstart=3000;
            TrialParameters=fminsearch(@(lambda)(fitblackbody(lambda,xx,yy)),bbstart,options);
        case 48
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewexpulse(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        otherwise
    end % switch peakshape

    % Check variables
    % sizeNewstart=size(newstart) % uncomment for testing
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
    switch peakshape(1)
        case 1
            A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 2
            A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 3
            A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 4
            A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 5
            A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 6
            A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 7
            A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 8
            A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
        case 9
            A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 10
            A(m,:)=upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 11
            A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
        case 13
            A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 14
            A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 15
            A(m,:)=BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 16
            A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 17
            A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 18
            A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 19
            A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 20
            A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 21
            A(m,:)=triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 22
            A(m,:)=peakfunction(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m));        
        case 23
            A(m,:)=downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));        
        case 24
            A(m,:)=nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 25
            A(m,:)=lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 26
            A(m,:)=linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 27
            A(m,:)=d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 28
            A(m,:)=polynomial(xx,coeff);
        case 29
            A(m,:)=segmented(xx,yy,PEAKHEIGHTS);
        case 30
            A(m,:)=voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 31
            A(m,:)=expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),-TrialParameters(3*m));        
        case 32
            A(m,:)=pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 33
            A(m,:)=GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
        case 34
             width(m)=abs(FIXEDPARAMETERS(m));
%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2))
            A(m,:)=voigt(xx,TrialParameters(m), width(m),extra);
        case 35
            A(m,:)=GL(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
        case 36
            A(m,:)=expgaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
        case 37
            A(m,:)=pearson(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
        case 38
            A(m,:)=explorentzian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),-TrialParameters(3*m));        
        case 40
            A(m,:)=sine(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 41
            A(m,:)=rectangle(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 42
            A(m,:)=ngaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 43
            A(m,:)=Gompertz(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 44
            A(m,:)=OneMinusExp(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 45
            A(m,:)=FourPL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 46
            A(m,:)=quadslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 47
            A(m,:)=blackbody(xx,TrialParameters(m));
        case 48
            A(m,:)=exppulse(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
    end % switch
    xxrange=max(xx)-min(xx);
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)+(xxrange.*(delta*randn)./(NumPeaks+1));
        newstart(parameter+1)=newstart(parameter+1)*(1+delta*(rand-.5)/100);
    end
end % for NumPeaks
% newstart=newstart; % <<<<<<<<<< % error check
% Multiplies each row by the corresponding amplitude and adds them up
if peakshape(1)==29, % Segmented linear
    model=segmented(xx,yy,PEAKHEIGHTS);
    TrialParameters=PEAKHEIGHTS;
    Heights=ones(size(PEAKHEIGHTS));
else
    if AUTOZERO==3,
        baseline=PEAKHEIGHTS(1);
        Heights=PEAKHEIGHTS(2:1+NumPeaks);
        model=Heights'*A+baseline;
    else
%          sizePeakHeights=size(PEAKHEIGHTS) %  % uncomment for testing
%          SizeA=size(A)  % uncomment for testing
        model=PEAKHEIGHTS'*A;
        Heights=PEAKHEIGHTS;
        baseline=0;
    end
end
if peakshape(1)==28, % polynomial;
    model=polynomial(xx,coeff);
    TrialParameters=PEAKHEIGHTS;
    Heights=ones(size(PEAKHEIGHTS));
end
% Compare trial model to data segment and compute the fit error
    MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=-BIPOLAR*10^100,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
%  ErrorVector(k)=MeanFitError; %  % uncomment for testing
end % for k (NumTrials)
    Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)));
    SStot=sum((yy-mean(yy)).^2);
    SSres=sum((yy-bestmodel).^2);
    Rsquared=1-(SSres./SStot);
    GOF=[LowestError Rsquared];
% Uncomment following 4 lines to monitor trail fit starts and errors.
% StartMatrix=StartMatrix;
% ErrorVector=ErrorVector;
% matrix=[StartMatrix ErrorVector']
% std(StartMatrix)
% Construct model from best-fit parameters
AA=zeros(NumPeaks,600);
xxx=linspace(min(xx),max(xx),600);
% xxx=linspace(min(xx)-length(xx),max(xx)+length(xx),200);
for m=1:NumPeaks,
   switch peakshape(1)
    case 1
        AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 6
        AA(m,:)=gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 7
        AA(m,:)=lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 8
        AA(m,:)=expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra*length(xxx)./length(xx))';
    case 9
        AA(m,:)=exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 10
        AA(m,:)=upsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 11
        AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDPARAMETERS(m));
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDPARAMETERS(m));
    case 13
        AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 14
        AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 15
        AA(m,:)=BWF(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 16
        AA(m,:)=gaussian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 17
        AA(m,:)=lorentzian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 18
        AA(m,:)=explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 19
        AA(m,:)=alphafunction(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 20
        AA(m,:)=voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 21
        AA(m,:)=triangular(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 22
        AA(m,:)=peakfunction(shapesvector(m),xxx,FitParameters(2*m-1),FitParameters(2*m),extra(m));        
    case 23
        AA(m,:)=downsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 24
        AA(m,:)=nbinpdf(xxx,FitParameters(2*m-1),FitParameters(2*m));    
    case 25
        AA(m,:)=lognormal(xxx,FitParameters(2*m-1),FitParameters(2*m));    
    case 26
        AA(m,:)=linslope(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 27
        AA(m,:)=d1gauss(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 28
        AA(m,:)=polynomial(xxx,coeff);
    case 29
    case 30
        AA(m,:)=voigt(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));        
    case 31
        AA(m,:)=expgaussian(xxx,FitParameters(3*m-2),FitParameters(3*m-1),-FitParameters(3*m)*length(xxx)./length(xx));        
    case 32
        AA(m,:)=pearson(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));        
    case 33
        AA(m,:)=GL(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m)); 
    case 34
                  width(m)=abs(FIXEDPARAMETERS(m));
%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 +
%                 gD(m).^2));
        AA(m,:)=voigt(xxx,FitParameters(m),width(m),extra);
    case 35
        AA(m,:)=GL(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);
    case 36
        AA(m,:)=expgaussian(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);    
    case 37
        AA(m,:)=pearson(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);
    case 38
        AA(m,:)=explorentzian(xxx,FitParameters(3*m-2),FitParameters(3*m-1),-FitParameters(3*m)*length(xxx)./length(xx));        
    case 40
        AA(m,:)=sine(xx,FitParameters(2*m-1),FitParameters(2*m));
    case 41
        AA(m,:)=rectangle(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 42
        AA(m,:)=ngaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 43
        AA(m,:)=Gompertz(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));  
    case 44
        AA(m,:)=OneMinusExp(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 45
        AA(m,:)=FourPL(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));    
    case 46
        AA(m,:)=quadslope(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 47
        AA(m,:)=blackbody(xxx,FitParameters(m));
    case 48
        AA(m,:)=exppulse(xxx,FitParameters(m),FitParameters(NumPeaks+1));
       otherwise
  end % switch
end % for NumPeaks

% Multiplies each row by the corresponding amplitude and adds them up
if peakshape(1)==29, % Segmented linear
    mmodel=segmented(xx,yy,PEAKHEIGHTS);
    baseline=0;
else
    heightsize=size(height');
    AAsize=size(AA);
    if heightsize(2)==AAsize(1),
        mmodel=height'*AA+baseline;
    else
        mmodel=height*AA+baseline;
    end
end
% Top half of the figure shows original signal and the fitted model.
if plots,
    subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
    hold on
end
if peakshape(1)==28, % Polynomial
     yi=polynomial(xxx,coeff);
else
    for m=1:NumPeaks,
        if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),end  % Plot the individual component peaks in green lines
        area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
        yi(m,:)=height(m)*AA(m,:); % Place y values of individual model peaks into matrix yi
    end
end
xi=xxx+xoffset; % Place the x-values of the individual model peaks into xi
residual=yy-bestmodel;

if plots,
    % Mark starting peak positions with vertical dashed magenta lines
    if peakshape(1)==16||peakshape(1)==17
    else
        if peakshape(1)==29, % Segmented linear
%            subplot(2,1,1);plot([PEAKHEIGHTS' PEAKHEIGHTS'],[0 max(yy)],'m--')
        else
            for marker=1:NumPeaks,
                markx=BestStart((2*marker)-1);
%                subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
            end % for
        end
    end % if peakshape

    % Plot the total model (sum of component peaks) in red lines
    if peakshape(1)==29, % Segmented linear
        mmodel=segmented(xx,yy,PEAKHEIGHTS);
       plot(xx+xoffset,mmodel,'r');  
    else
       plot(xxx+xoffset,mmodel,'r');  
    end
    hold off;
    lyy=min(yy);
    uyy=max(yy)+(max(yy)-min(yy))/10;
    if BIPOLAR,
        axis([min(xx) max(xx) lyy uyy]);
        ylabel('+ - mode')
    else
        axis([min(xx) max(xx) 0 uyy]);
        ylabel('+ mode')
    end
    switch AUTOZERO,
        case 0
            title(['peakfit.m Version 8   No baseline correction'])
        case 1
            title(['peakfit.m Version 8   Linear baseline subtraction'])
        case 2
            title(['peakfit.m Version 8   Quadratic subtraction baseline'])
        case 3
            title(['peakfit.m Version 8   Flat baseline correction'])
    end
 
    switch peakshape(1)
        case {4,20,34,37,42}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Shape Constant = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
        case {5,8,18,36}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Time Constant = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
        case 13
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      % Gaussian = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
        case {14,15,22,35}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      extra = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
        case 28
            xlabel(['Shape = ' ShapeString '      Order = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(1000*LowestError)/1000) ] )
        case 43
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '        Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
         otherwise
            if peakshape(1)==29, % Segmented linear
                xlabel(['Breakpoints = ' num2str(NumPeaks) '     Shape = ' ShapeString  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
            else
                xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
            end % if peakshape(1)==29
    end % switch peakshape(1)

    % Bottom half of the figure shows the residuals and displays RMS error
    % between original signal and model
    % residual=yy-bestmodel;
    subplot(2,1,2);plot(xx+xoffset,residual,'m.')
    axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
    xlabel('Residual Plot')
    if NumTrials>1,
       title(['Best of ' num2str(NumTrials) ' fits'])
    else
       title(['Single fit'])
    end
end % if plots

% Put results into a matrix FitResults, one row for each peak, showing peak index number,
% position, amplitude, and width.
FitResults=zeros(NumPeaks,6);
%  FitParameters=FitParameters
switch peakshape(1),
    case {6,7,8,48}, % equal-width peak models only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            end
        end
    case {11,12,34,35,36,37}, % Fixed-width shapes only
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(m));
            if peakshape==34,
                gD(m)=width(m);
                gL(m)=extra.*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            end
        end
    case {16,17}, % Fixed-position shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
            end
        end
    case 28,   % Simple polynomial fit
        FitResults=PEAKHEIGHTS;
    case 29, % Segmented linear fit
        FitResults=PEAKHEIGHTS;
    case {30,31,32,33,38,43} % Special case of shapes with 3 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(3*m-1));
            if peakshape==30,
                gD(m)=width(m);
                gL(m)=FitParameters(3*m).*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) width(m) area(m)] FitParameters(3*m)];
            end
        end
    case 47 % Shapes with 1 iterated variable
        FitParameters=FitParameters;
    otherwise % Normal shapes with 2 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(2*m));
            if peakshape==20,
                gD=width(m);
                gL=extra.*gD;
                width(m) = 2.*(0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)]];
            end % if m==1
            
        end % for m=1:NumPeaks,
end % switch peakshape(1)
  
% Rearrange fit results for Gompertz to Bo, Kh, and L
if peakshape(1)==43;
    for m=1:NumPeaks,
        FitResults(m,2)=FitResults(m,2).*FitResults(m,3);
        FitResults(m,4)=FitResults(m,3).*FitResults(m,4);
        FitResults(m,3)=1;
    end
end

%Sort FitResults
FitResults=sortrows(FitResults,2);

% Display Fit Results on lower graph
if plots,
    % Display Fit Results on lower  graph
    subplot(2,1,2);
    startx=min(xx)+(max(xx)-min(xx))./20;
    dxx=(max(xx)-min(xx))./10;
    dyy=((max(residual)-min(residual))./10);
    starty=max(residual)-dyy;
    FigureSize=get(gcf,'Position');
    switch peakshape(1)
        case {9,19,10,23,40}  % Pulse and sigmoid shapes only
            text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
        case 28, % Polynomial
            text(startx,starty+dyy/2,['Polynomial coefficients'] );
        case 29 % Segmented linear
             text(startx,starty+dyy/2,['x-axis breakpoints'] );
        case {30,31,32,33,38} % Special case of shapes with 3 iterated variables
            text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area       Shape factor'] );            
        case 43 % 3 parameter Gompertz
            text(startx,starty+dyy/2,['Peak #           Bo             Height            Kh                Area                 L'] );            
        otherwise
            text(startx,starty+dyy/2,['Peak #          Position         Height         Width             Area   '] );
    end
    % Display FitResults using sprintf
    if peakshape(1)==28||peakshape(1)==29, % Polynomial or segmented linear
        for number=1:length(FitResults),
            column=1;
            itemstring=sprintf('%0.4g',FitResults(number));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-number.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,['                ' itemstring]);
        end
    else
        for peaknumber=1:NumPeaks,
            for column=1:5,
                itemstring=sprintf('%0.4g',FitResults(peaknumber,column));
                xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
                yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
                text(xposition,yposition,itemstring);
            end
        end
        xposition=startx;
        yposition=starty-(peaknumber+1).*dyy.*(400./FigureSize(4));
        if AUTOZERO==3,
            text(xposition,yposition,[ 'Baseline= ' num2str(baseline) ]);
        end % if AUTOZERO
    end % if peakshape(1)
    if peakshape(1)==30 || peakshape(1)==31 || peakshape(1)==32 || peakshape(1)==33 || peakshape(1)==38 || peakshape(1)==43,
        for peaknumber=1:NumPeaks,
            column=6;
            itemstring=sprintf('%0.4g',FitParameters(3*peaknumber));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,itemstring);
        end
    end
end % if plots

if NumArgOut==8,
    if plots,disp('Computing bootstrap sampling statistics.....'),end
    BootstrapResultsMatrix=zeros(6,100,NumPeaks);
    BootstrapErrorMatrix=zeros(1,100,NumPeaks);
    clear bx by
    tic;
    for trial=1:100,
        n=1;
        bx=xx;
        by=yy;
        while n<length(xx)-1,
            if rand>.5,
                bx(n)=xx(n+1);
                by(n)=yy(n+1);
            end
            n=n+1;
        end
        bx=bx+xoffset;
        [FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDPARAMETERS,shapesvector);
        for peak=1:NumPeaks,
            switch peakshape(1)
                case {30,31,32,33,38,43}
                    BootstrapResultsMatrix(1:6,trial,peak)=FitResults(peak,1:6);
                otherwise
                    BootstrapResultsMatrix(1:5,trial,peak)=FitResults(peak,1:5);
            end
            BootstrapErrorMatrix(:,trial,peak)=BootFitError;
        end
    end
    if plots,toc;end
    for peak=1:NumPeaks,
        if plots,
            disp(' ')
            % Label columns for bootstrap results
            switch peakshape(1)
                case {9,19,10,23,40}  % Pulse and sigmoid shapes only
                    disp(['Peak #',num2str(peak) '         tau1           Height         tau2           Area'] );
                case {30,31,32,33,38} % Special case of shapes with 3 iterated variables
                    disp(['Peak #',num2str(peak) '         Position        Height         Width             Area       Shape factor'] );
                case 43 % 3 parameter Gompertz
                    disp(['Peak #',num2str(peak) '         Bo             Height         Kh                       L'] );
                otherwise
                    disp(['Peak #',num2str(peak) '         Position         Height         Width             Area'] );
            end
            
        end % if plots
        BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
        BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
        BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
        PercentRSD=100.*BootstrapSTD./BootstrapMean;
        PercentIQR=100.*BootstrapIQR./BootstrapMean;
        BootstrapMean=BootstrapMean(2:6);
        BootstrapSTD=BootstrapSTD(2:6);
        BootstrapIQR=BootstrapIQR(2:6);
        PercentRSD=PercentRSD(2:6);
        PercentIQR=PercentIQR(2:6);
        format short g
        if plots,
            disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
            disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
            disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
            disp(['Percent RSD:    ', num2str(PercentRSD)])
            disp(['Percent IQR:    ', num2str(PercentIQR)])
        end % if plots
        BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR];
    end % peak=1:NumPeaks,
end % if NumArgOut==8,
if AUTOZERO==3;
else
    baseline=bkgcoef;
end
% ----------------------------------------------------------------------
function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters,shapesvector)
% Based on peakfit Version 3: June, 2012. 
global PEAKHEIGHTS FIXEDPARAMETERS BIPOLAR MINWIDTH coeff
format short g
format compact
warning off all
FIXEDPARAMETERS=fixedparameters;
xoffset=0;
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;
coeff=0;
LOGPLOT=0;

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.000001,'Display','off','MaxFunEvals',1000 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials,
    % StartVector=newstart
    switch peakshape(1)
        case 1
            TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 2
            TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 3
            TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 4
            TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 5
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 6
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 7
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewlorentzian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 8
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 9
            TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 10
            TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 11
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 12
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 13
            TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 14
            TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 15
            TrialParameters=fminsearch(@(lambda)(fitBWF(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 16
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 17
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 18
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
        case 19
            TrialParameters=fminsearch(@(lambda)(alphafunction(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 20
            TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 21
            TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 22
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,NumPeaks,shapesvector,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 23
            TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstart,optionst);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 24
            TrialParameters=fminsearch(@(lambda)(fitnbinpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 25
            TrialParameters=fminsearch(@(lambda)(fitlognpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 26
            TrialParameters=fminsearch(@(lambda)(fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options);
        coeff=TrialParameters;
        case 27
            TrialParameters=fminsearch(@(lambda)(fitd1gauss(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 28
            TrialParameters=fitpolynomial(xx,yy,extra);
        case 29
            TrialParameters=fminsearch(@(lambda)(fitsegmented(lambda,xx,yy)),newstart,options);
        case 30
            TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),newstart);
        case 31
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussianv(lambda,zxx,zyy)),newstart);
        case 32
            TrialParameters=fminsearch(@(lambda)(fitpearsonv(lambda,xx,yy)),newstart);
        case 33
            TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart);
        case 34
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWVoigt(lambda,xx,yy,extra)),fixedstart,options);
        case 35
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWGL(lambda,xx,yy,extra)),fixedstart,options);
        case 36
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWExpGaussian(lambda,xx,yy,extra)),fixedstart,options);
        case 37
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWPearson(lambda,xx,yy,extra)),fixedstart,options);
        case 38
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzianv(lambda,zxx,zyy)),newstart);         
        case 40
             TrialParameters=fminsearch(@(lambda)(fitsine(lambda,xx,yy)),newstart,options);
        case 41
            TrialParameters=fminsearch(@(lambda)(fitrectangle(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 42
            TrialParameters=fminsearch(@(lambda)(fitngaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 43
            TrialParameters=fminsearch(@(lambda)(fitGompertz(lambda,xx,yy)),newstart);
        case 44
            TrialParameters=fminsearch(@(lambda)(fitOneMinusExp(lambda,xx,yy)),newstart);           
        case 45
            TrialParameters=fminsearch(@(lambda)(fitFourPL(lambda,xx,yy)),newstart);  
        case 46
            TrialParameters=fminsearch(@(lambda)(fitquadslope(lambda,xx,yy)),newstart,options);
        case 47
            bbstart=3000;
            TrialParameters=fminsearch(@(lambda)(fitblackbody(lambda,xx,yy)),bbstart,options);
        case 48
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewexpulse(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        otherwise
    end % switch peakshape
    
for peaks=1:NumPeaks,
     peakindex=2*peaks-1;
     newstart(peakindex)=start(peakindex)-xoffset;
end

    % Construct model from Trial parameters
    A=zeros(NumPeaks,n);
    for m=1:NumPeaks,
        switch peakshape(1)
            case 1
                A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 2
                A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 3
                A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 4
                A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 5
                A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
            case 6
                A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
            case 7
                A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
            case 8
                A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
            case 9
                A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 10
                A(m,:)=upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 11
                A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
            case 12
                A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
            case 13
                A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 14
                A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 15
                A(m,:)=BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 16
                A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
            case 17
                A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
            case 18
                A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
            case 19
                A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 20
                A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 21
                A(m,:)=triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 22
                A(m,:)=peakfunction(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m));
            case 23
                A(m,:)=downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));      
            case 24
                A(m,:)=nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 25
                A(m,:)=lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 26
                A(m,:)=linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 27
                A(m,:)=d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m));       
            case 28
                A(m,:)=polynomial(xx,TrialParameters(2*m-1),TrialParameters(2*m));       
            case 29
                A(m,:)=segmented(xx,yy,PEAKHEIGHTS);
            case 30
                A(m,:)=voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 31
                A(m,:)=expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 32
                A(m,:)=pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 33
                A(m,:)=GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
            case 34
                width(m)=abs(FIXEDPARAMETERS(m));

%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2))

                A(m,:)=voigt(xx,TrialParameters(m), width(m),extra);
            case 35
                A(m,:)=GL(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
            case 36
                A(m,:)=expgaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
            case 37
                A(m,:)=pearson(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
            case 38
                A(m,:)=explorentzian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));                    
            case 40
                A(m,:)=sine(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 41
                A(m,:)=rectangle(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 42
                A(m,:)=ngaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 43
                A(m,:)=Gompertz(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
            case 44
                A(m,:)=OneMinusExp(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 45
                A(m,:)=FourPL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 46
                A(m,:)=quadslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 47
                A(m,:)=blackbody(xx,TrialParameters(m));
            case 48
                A(m,:)=exppulse(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        end % switch
    end % for
    
    % Multiplies each row by the corresponding amplitude and adds them up
    if peakshape(1)==29, % Segmented linear
        model=segmented(xx,yy,PEAKHEIGHTS);
        TrialParameters=coeff;
        Heights=ones(size(coeff));
    else
        if AUTOZERO==3,
            baseline=PEAKHEIGHTS(1);
            Heights=PEAKHEIGHTS(2:1+NumPeaks);
            model=Heights'*A+baseline;
        else
            model=PEAKHEIGHTS'*A;
            Heights=PEAKHEIGHTS;
            baseline=0;
        end
    end
    
    % Compare trial model to data segment and compute the fit error
    MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
    % Take only the single fit that has the lowest MeanFitError
    if MeanFitError<LowestError,
        if min(Heights)>=-BIPOLAR*10^100,  % Consider only fits with positive peak heights
            LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
            FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
            height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        end % if min(PEAKHEIGHTS)>0
    end % if MeanFitError<LowestError
end % for k (NumTrials)
    Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)));
    SStot=sum((yy-mean(yy)).^2);
    SSres=sum((yy-bestmodel).^2);
    Rsquared=1-(SSres./SStot);
    GOF=[LowestError Rsquared];
for m=1:NumPeaks,
    area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end

% Put results into a matrix FitResults, one row for each peak, showing peak index number,
% position, amplitude, and width.
FitResults=zeros(NumPeaks,6);
switch peakshape(1),
    case {6,7,8}, % equal-width peak models only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            end
        end
    case {11,12,35,36,37}, % Fixed-width shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            end
        end

    case {16,17}, % Fixed-position shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
            end
        end
    case 28,   % Simple polynomial fit
        FitResults=PEAKHEIGHTS;
    case 29, % Segmented linear fit
        FitResults=PEAKHEIGHTS;
    case {30,31,32,33,38,43} % Special case of shapes with 3 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(3*m-1));
            if peakshape==30,
                gD(m)=width(m);
                gL(m)=FitParameters(3*m).*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)]];
            end
        end
    otherwise % Normal shapes with 2 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(2*m));
            if peakshape==20,
                gD=width(m);
                gL=extra.*gD;
                width(m) = 2.*(0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)]];
            end % if m==1
        end % for m=1:NumPeaks,
end % switch peakshape(1)
if peakshape==34,
        DW=2*(0.5346*a*1.2772 + sqrt(0.2166*a*1.2772.^2 + 1.2772.^2))
end
% ----------------------------------------------------------------------
function start=calcstart(xx,NumPeaks,xoffset)
  n=max(xx)-min(xx);
  start=[];
  startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
  for marker=1:NumPeaks,
      markx=startpos(marker)+ xoffset;
      start=[start markx n/ (3.*NumPeaks)];
  end % for marker
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function err = fitgaussian(lambda,t,y)
% Fitting function for multiple Gaussian peaks.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
%    if lambda(2*j)<MINWIDTH,lambda(2*j)=MINWIDTH;end
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitngaussian(lambda,t,y,shapeconstant)
%   Fitting functions for multiple flattened Gaussian peaks.
% T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = ngaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitewgaussian(lambda,t,y)
% Fitting function for multiple Gaussian peaks with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFWGaussian(lambda,t,y)
%	Fitting function for multiple fixed width Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),FIXEDPARAMETERS(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFPGaussian(lambda,t,y)
%	Fitting function for multiple fixed-position Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,FIXEDPARAMETERS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFPLorentzian(lambda,t,y)
%	Fitting function for multiple fixed-position Lorentzians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,FIXEDPARAMETERS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWLorentzian(lambda,t,y)
%	Fitting function for multiple fixed width Lorentzians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),FIXEDPARAMETERS(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitewlorentzian(lambda,t,y)
% Fitting function for multiple  Lorentzian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.60056120439323.*wid)).^2);
% ----------------------------------------------------------------------
function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = flattened Gaussian centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
% Shape is Gaussian when n=1. Becomes more rectangular as n increases.
%  T. C. O'Haver, 1988, revised 2014
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
if n>0,
    g = 1-(10.^-(n.*gaussian(x,pos,wid)));
    g=g./max(g);
else
    g = gaussian(x,pos,wid);
end
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for multiple lorentzians, lambda(1)=position, lambda(2)=width
%	Fitgauss assumes a lorentzian function 
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% [lambda PEAKHEIGHTS err]
% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
% ----------------------------------------------------------------------
function err = fitlogistic(lambda,t,y)
%	Fitting function for multiple logistic peaks, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlogistic assumes a logistic function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991 
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function err = fittriangular(lambda,t,y)
%	Fitting function for multiple triangular, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fittriangular assumes a triangular function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = triangular(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = triangular(x,pos,wid)
%Triangle function.  pos=position; wid=half-width (both scalar)
%trianglar(x,pos,wid), where x may be scalar or vector,
%pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
% Example
% x=[0:.1:10];plot(x,trianglar(x,5.5,2.3),'.')
g=1-(1./wid) .*abs(x-pos);
for i=1:length(x),  
if g(i)<0,g(i)=0;end
end
% ----------------------------------------------------------------------
function err = fitrectangle(lambda,t,y)
%	Fitting function for multiple rectangle, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitrectangle assumes a rectangle function 
%  T. C. O'Haver, May 2016
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = rectangle(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = rectangle(x,pos,wid)
%rectangle function.  pos=position; wid=half-width (both scalar)
%rectangle(x,pos,wid), where x may be scalar or vector,
%pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 2016
% Example
% x=[0:.1:10];plot(x,rectangle(x,5.5,2.3),'.')
g=zeros(size(x));
hw=wid./2;
for i=1:length(x),  
if x(i)<pos-hw,g(i)=0;end
if x(i)>pos-hw,g(i)=1;end
if x(i)>pos+hw,g(i)=0;end
end
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for multiple Pearson 7 bands.
% T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitpearsonv(lambda,t,y)
% Fitting functions for multiple pearson functions with independently variable
% percent Gaussian
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = pearson(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitFWPearson(lambda,t,y,shapeconstant)
%	Fitting function for multiple fixed width Pearson7
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = pearson(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function. 
% g = pearson(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
%  T. C. O'Haver, 1990  
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function err = fitexpgaussian(lambda,t,y,timeconstant)
%   Fitting functions for multiple exponentially-broadened Gaussian bands signal.
%  T. C. O'Haver, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexplorentzian(lambda,t,y,timeconstant)
%   Fitting functions for multiple exponentially-broadened lorentzian band signal.
%  T. C. O'Haver, 2013.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexplorentzianv(lambda,t,y)
% Fitting functions for multiple exponentially-broadened Gaussians with
% independently variable time constants
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = explorentzian(t,lambda(3*j-2),lambda(3*j-1),-lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexpewgaussian(lambda,t,y,timeconstant)
% Fitting function for multiple exponentially-broadened Gaussian bands with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexpgaussianv(lambda,t,y)
% Fitting functions for multiple exponentially-broadened Gaussians with
% independently variable time constants
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = expgaussian(t,lambda(3*j-2),lambda(3*j-1),-lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitFWExpGaussian(lambda,t,y,shapeconstant)
%	Fitting function for multiple fixed width exponentially-broadened gaussian
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2006
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function g = explorentzian(x,pos,wid,timeconstant)
%  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2013
g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant t by multiplying Fourier transforms and inverse
% transforming the result.
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
% ----------------------------------------------------------------------
function err = fitexppulse(tau,x,y)
% Iterative fit of the sum of exponential pulses
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitewexpulse(lambda,t,y)
% Fitting function for multiple Gaussian peaks with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = exppulse(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form 
% g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* (p>0);
g = p';
% ----------------------------------------------------------------------
function err = fitalphafunction(tau,x,y)
% Iterative fit of the sum of alpha funciton
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = alphafunction(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = alphafunction(x,pos,spoint)
% alpha function.  pos=position; wid=half-width (both scalar)
% alphafunction(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% Taekyung Kwon, July 2013  
g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
for m=1:length(x);if g(m)<0;g(m)=0;end;end
% ----------------------------------------------------------------------
function err = fitdownsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of multiple
% downward moving sigmiods 
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = downsigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitupsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of multiple
% upwards moving sigmiods
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = upsigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=downsigmoid(x,t1,t2)
 % down step sigmoid
g=.5-.5*erf(real((x-t1)/sqrt(2*t2)));
% ----------------------------------------------------------------------
function g=upsigmoid(x,t1,t2)
% up step sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); 
% ----------------------------------------------------------------------
function err = fitGL(lambda,t,y,shapeconstant)
%   Fitting functions for multiple Gaussian/Lorentzian blend peaks
% T. C. O'Haver (toh@umd.edu), 2012.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitFWGL(lambda,t,y,shapeconstant)
%	Fitting function for a multiple fixed width Gaussian/Lorentzian blend
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = GL(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitGLv(lambda,t,y)
% Fitting functions for multiple Gaussian/Lorentzian blend functions with
% independently variable percent Gaussian
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = GL(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend. m = percent Gaussian character
% pos=position; wid=half-width
% m = percent Gaussian character.
%  T. C. O'Haver, 2012
% sizex=size(x)
% sizepos=size(pos)
% sizewid=size(wid)
% sizem=size(m)
g=2.*((m/100).*gaussian(x,pos,wid)+(1-(m(1)/100)).*lorentzian(x,pos,wid))/2;
% ----------------------------------------------------------------------
function err = fitvoigt(lambda,t,y,shapeconstant)
% Fitting functions for multiple Voigt profile function
% T. C. O'Haver (toh@umd.edu), 2013.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitFWVoigt(lambda,t,y,shapeconstant)
%	Fitting function for multiple fixed width Voigt
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
% numpeaksfitFWVoigt=numpeaks
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = voigt(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitvoigtv(lambda,t,y)
% Fitting functions for multiple Voigt profile function with independently variable
% alphas
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = voigt(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=voigt(xx,pos,gD,alpha)
% Voigt profile function. xx is the independent variable (energy,
% wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
% shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
% Based on Chong Tao's "Voigt lineshape spectrum simulation", 
% File ID: #26707
% alpha=alpha
gL=alpha.*gD;
gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2);
x = gL/gV;
% sizeabs=size(abs(xx-pos))
% sizegV=size(gV)
y = abs(xx-pos)./gV;
g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
g=g./max(g);
% ----------------------------------------------------------------------
function err = fitBiGaussian(lambda,t,y,shapeconstant)
%   Fitting functions for multiple BiGaussian.
% T. C. O'Haver (toh@umd.edu),  2012.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = BiGaussian(x,pos,wid,m)
% BiGaussian (different width on leading edge and trailing edge).
% pos=position; wid=width 
% m = ratio of widths of right-hand to left-hand halves.
% If m=1, it becomes identical to a Gaussian.
% Verison 2, T. C. O'Haver, 2012
% Example: Plots Gaussian and BiGaussian with m=3, both with halfwidth=20.
% x=[1:100];
% g=gaussian(x,50,20);
% bg=bigaussian(x,50,20,3);
% plot(x,g,x,bg)
%
lx=length(x);
hx=val2ind(x,pos);
hwid=2.*wid;
g(1:hx)=gaussian(x(1:hx),pos,hwid/(m+1));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,m.*hwid/(m+1));
% ----------------------------------------------------------------------
function err = fitBWF(lambda,t,y,shapeconstant)
%   Fitting function for multiple Breit-Wigner-Fano.
% T. C. O'Haver (toh@umd.edu),  2014.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BWF(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = BWF(x,pos,wid,m)
% BWF (Breit-Wigner-Fano) http://en.wikipedia.org/wiki/Fano_resonance
% pos=position; wid=width; m=Fano factor
%  T. C. O'Haver, 2014
y=((m*wid/2+x-pos).^2)./(((wid/2).^2)+(x-pos).^2);
% y=((1+(x-pos./(m.*wid))).^2)./(1+((x-pos)./wid).^2);
g=y./max(y);
% ----------------------------------------------------------------------
function err = fitnbinpdf(tau,x,y)
% Fitting function for iterative fit to the sum of multiple
% Negative Binomial Distributions
% (http://www.mathworks.com/help/stats/negative-binomial-distribution.html)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = nbinpdf(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitlognpdf(tau,x,y)
% Fitting function for iterative fit to the sum of multiple
% Lognormal Distributions
% (http://www.mathworks.com/help/stats/lognormal-distribution.html)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = lognormal(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991  
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitsine(tau,x,y)
% Fitting function for iterative fit to the sum of multiple
% sine waves (alpha test, NRFPT)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = sine(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=sine(t,f,phase) 
% Sine wave (alpha test)
g=sin(2*pi*f*(t+phase));
% ----------------------------------------------------------------------
function err = fitd1gauss(lambda,t,y)
%   Fitting functions for multiple first derivative of a Gaussian
%  T. C. O'Haver, 2014
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = d1gauss(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function y=d1gauss(x,p,w)
% First derivative of Gaussian (alpha test)
y=-(5.54518.*(x-p).*exp(-(2.77259.*(p-x).^2)./w^2))./w.^2;
y=y./max(y);
% ----------------------------------------------------------------------
function coeff = fitpolynomial(t,y,order)
coeff=polyfit(t,y,order);
% order=order
% coeff=coeff
% ----------------------------------------------------------------------
function y=polynomial(t,coeff)
y=polyval(coeff,t);
% ----------------------------------------------------------------------
function err = fitsegmented(lambda,t,y)
%   Fitting functions for articulated segmented linear
%  T. C. O'Haver, 2014
global LOGPLOT
breakpoints=[t(1) lambda max(t)];
z = segmented(t,y,breakpoints);
% lengthz=length(z);
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y);
end
% ----------------------------------------------------------------------
function yi=segmented(x,y,segs)
global PEAKHEIGHTS
clear yy
for n=1:length(segs)
  yind=val2ind(x,segs(n));
  yy(n)=y(yind(1));
end
yi=interp1(segs,yy,x);
PEAKHEIGHTS=segs;
% ----------------------------------------------------------------------
function err = fitlinslope(tau,x,y)
% Fitting function for iterative fit to linear function
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    z = (x.*tau(2*j-1)+tau(2*j))';
    A(:,j) = z./max(z);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=linslope(x,slope,intercept)
y=x.*slope+intercept;
% y=y./max(y);
% ----------------------------------------------------------------------
function err = fitquadslope(lambda,x,y)
% Fitting function for iterative fit to linear function
global PEAKHEIGHTS LOGPLOT
A = zeros(length(x),round(length(lambda)/2));
for j = 1:length(lambda)/2,
     A(:,j)=quadslope(x,lambda(1),lambda(2));
end
PEAKHEIGHTS=A\y';
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=quadslope(x,boa,coa) % normalized quadratic
 y=(x.^2+(boa).*x+coa);
y=y./max(y);
% ----------------------------------------------------------------------
function err = fitGompertz(lambda,t,y)
% Fitting functions for multiple Gompertz function
% T. C. O'Haver (toh@umd.edu), 2016.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
% Sizelambda=size(lambda)
for j = 1:length(lambda)/3,
    A(:,j) = Gompertz(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
% PEAKHEIGHTS=1;
% SizeA=size(A)
% sizwPEAKHEIGHTS=size(PEAKHEIGHTS)
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=Gompertz(t,Bo,Kh,L)
% A Gompertz curve or Gompertz function, named after Benjamin Gompertz, is
% a sigmoid function. It is a type of mathematical model for a time series,
% where growth is slowest at the start and end of a time period. The
% right-hand or future value asymptote of the function is approached much
% more gradually by the curve than the left-hand or lower valued asymptote,
% in contrast to the simple logistic function in which both asymptotes are
% approached by the curve symmetrically. It is a special case of the
% generalized logistic function.
% 
% Example:
% x=1:.1:10;y=gompertz(x,6,3,4);plot(x,y)
y=Bo*exp(-exp((Kh*exp(1)/Bo)*(L-t) +1));
% ----------------------------------------------------------------------
function err = fitFourPL(lambda,t,y)
% Fitting functions for multiple Gompertz function
% T. C. O'Haver (toh@umd.edu), 2016.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = FourPL(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
PEAKHEIGHTS=1;
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=FourPL(x,miny,slope,ip)
% Normalized four parameters logistic regression  
% https://psg.hitachi-solutions.com/masterplex/blog/the-4-parameter-logisti
% c-4pl-nonlinear-regression-model
% miny = minimum asymptote. In an ELISA assay where you have a standard
% curve, this can be thought of as the response value at 0 standard
% concentration. slope = Hill slope. The Hill Slope or slope factor refers
% to the steepness of the curve. It could either be positive or negative.
% As the absolute value of the Hill slope increases, so does the steepness
% of the curve. ip = inflection point: The inflection point is defined as the
% point on the curve where the curvature changes direction or signs. This
% can be better explained if you can imagine the concavity of a sigmoidal
% curve. The inflection point is where the curve changes from being concave
% upwards to concave downwards. maxy = maximum asymptote. In an ELISA assay
% where you have a standard curve, this can be thought of as the response
% value for infinite standard concentration.
% Example:
% x=0:20;
% miny=0;slope=5;ip=10;d=0;maxy=10;
% y=FourPL(x,miny,slope,ip,maxy);plot(x,y)
%
y = 1+(miny-1)./(1+(x./ip).^slope);
% ----------------------------------------------------------------------
function err = fitOneMinusExp(lambda,t,y)
%   Fitting functions for multiple first derivative of a Gaussian
%  T. C. O'Haver, 2014
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = OneMinusExp(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = OneMinusExp(x,pos,wid)
% OneMinusExp(x,pos,wid) = 1-exp(-wid.*(x-pos));
% Example:  x=0:10;y=1-exp(-(.5.*x));plot(x,y)
g = 1-exp(-wid.*(x-pos));
% ----------------------------------------------------------------------
function err = fitblackbody(lambda,wavelength,y)
%  Fitting function for a blackbody spectrum.
%  T. C. O'Haver, May 2008
global PEAKHEIGHTS AUTOZERO BIPOLAR
% sizelambda=size(lambda)
radiance = blackbody(wavelength,lambda(1));
% if AUTOZERO==3,A=[ones(size(y))' A];end
% if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
% sizePEAKHEIGHTS=size(PEAKHEIGHTS)
% sizey=size(y)
PEAKHEIGHTS = radiance/y;
z = radiance*PEAKHEIGHTS;
err = norm(z-y);
% ----------------------------------------------------------------------
function radiance=blackbody(wavelength,temperature)
radiance = 1.19111E+16*wavelength.^(-5)./(exp(14380000./(wavelength*temperature))-1);
% ----------------------------------------------------------------------
function b=iqr(a)
% b = IQR(a)  returns the interquartile range of the values in a.  For
%  vector input, b is the difference between the 75th and 25th percentiles
%  of a.  For matrix input, b is a row vector containing the interquartile
%  range of each column of a.
%  T. C. O'Haver, 2012
mina=min(a);
sizea=size(a);
NumCols=sizea(2);
for n=1:NumCols,b(:,n)=a(:,n)-mina(n);end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));
% ----------------------------------------------------------------------
function err = fitmultiple(lambda,xx,y,numpeaks,shapesvector,extra)
% Fitting function for a multiple-shape band signal.
% The sequence of peak shapes are defined by the vector "shape".
% The vector "extra" determines the shape of variable-shape peaks.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT coeff FIXEDPARAMETERS
% FIXEDPARAMETERS=FIXEDPARAMETERS % testing
% PEAKHEIGHTS=PEAKHEIGHTS % testing
% lambda=lambda
% numpeaks=round(length(lambda)/2)
% sizeshapesvector=size(shapesvector)
A=zeros(numpeaks,length(xx));
for j = 1:numpeaks,
    if shapesvector(j)==28,
        coeff=polyfit(xx,y,extra(j));
        A(j,:) = polyval(coeff,xx);
    else
        % sizeA=size(A)
        switch shapesvector(j)
            case 1
                A(j,:)=gaussian(xx,lambda(2*j-1),lambda(2*j));
            case 2
                A(j,:)=lorentzian(xx,lambda(2*j-1),lambda(2*j));
            case 3
                A(j,:)=logistic(xx,lambda(2*j-1),lambda(2*j));
            case 4
                A(j,:)=pearson(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 5
                A(j,:)=expgaussian(xx,lambda(2*j-1),lambda(2*j),-extra(j))';
            case 6
                A(j,:)=gaussian(xx,lambda(j),lambda(NumPeaks+1));
            case 7
                A(j,:)=lorentzian(xx,lambda(j),lambda(NumPeaks+1));
            case 8
                A(j,:)=expgaussian(xx,lambda(j),lambda(NumPeaks+1),-extra(j))';
            case 9
                A(j,:)=exppulse(xx,lambda(2*j-1),lambda(2*j));
            case 10
                A(j,:)=upsigmoid(xx,lambda(2*j-1),lambda(2*j));
            case 11
                A(j,:)=gaussian(xx,lambda(j),FIXEDPARAMETERS(j));
            case 12
                A(j,:)=lorentzian(xx,lambda(j),FIXEDPARAMETERS(j));
            case 13
                A(j,:)=GL(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 14
                A(j,:)=BiGaussian(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 15
                A(j,:)=BWF(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 16
                A(j,:)=gaussian(xx,FIXEDPOSITIONS(j),lambda(j));
            case 17
                A(j,:)=lorentzian(xx,FIXEDPOSITIONS(j),lambda(j));
            case 18
                A(j,:)=explorentzian(xx,lambda(2*j-1),lambda(2*j),-extra(j))';
            case 19
                A(j,:)=alphafunction(xx,lambda(2*j-1),lambda(2*j));
            case 20
                A(j,:)=voigt(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 21
                A(j,:)=triangular(xx,lambda(2*j-1),lambda(2*j));
            case 22
                A(j,:)=peakfunction(shapesvector(j),xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 23
                A(j,:)=downsigmoid(xx,lambda(2*j-1),lambda(2*j));      
            case 24
                A(j,:)=nbinpdf(xx,lambda(2*j-1),lambda(2*j));
            case 25
                A(j,:)=lognormal(xx,lambda(2*j-1),lambda(2*j));
            case 26
                A(j,:)=linslope(xx,lambda(2*j-1),lambda(2*j));
            case 27
                A(j,:)=d1gauss(xx,lambda(2*j-1),lambda(2*j));       
            case 28
                A(j,:)=polynomial(xx,lambda(2*j-1),lambda(2*j));       
            case 29
                A(j,:)=segmented(xx,yy,PEAKHEIGHTS);
            case 30
                A(j,:)=voigt(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
            case 31
                A(j,:)=expgaussian(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
            case 32
                A(j,:)=pearson(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
            case 33
                A(j,:)=GL(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
            case 34
                width(j)=abs(FIXEDPARAMETERS(j));
                 % [j lambda(j) width(j) extra(j)]
%                 gD(j)=width(j);
%                 gL(j)=extra.*gD(j);
%                 width(j) = 2.*(0.5346*gL(j) + sqrt(0.2166*gL(j).^2 +
%                 gD(j).^2))
                A(j,:)=voigt(xx,lambda(j),width(j),extra(j));
                % figure(5);plot(A(j,:));figure(1)
            case 35
                A(j,:)=GL(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
            case 36
                A(j,:)=expgaussian(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
            case 37
                A(j,:)=pearson(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
            case 38
                A(j,:)=explorentzian(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));                    
            case 40
                A(j,:)=sine(xx,lambda(2*j-1),lambda(2*j));
            case 41
                A(j,:)=rectangle(xx,lambda(2*    j-1),lambda(2*j));
            case 42
                A(j,:)=ngaussian(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 43
                A(j,:)=Gompertz(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
            case 44
                A(j,:)=OneMinusExp(xx,lambda(2*j-1),lambda(2*j));
            case 45
                A(j,:)=FourPL(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
            case 46
                A(j,:)=quadslope(xx,lambda(2*j-1),lambda(2*j));
            case 47
                A(j,:)=blackbody(xx,lambda(j));
            case 48
                A(j,:)=exppulse(xx,lambda(j),lambda(NumPeaks+1));
        end % switch
    end % if shapesvector
end % for j=1:numpeaks,
if AUTOZERO==3,A=[ones(size(y))' A];end
% sizeA=size(A)
% max(A')
% sizeyp=size(y')
if BIPOLAR,PEAKHEIGHTS=A'\y';else PEAKHEIGHTS=abs(A'\y');end
% PEAKHEIGHTS=PEAKHEIGHTS
z = A'*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function p=peakfunction(shape,x,pos,wid,extra,coeff)
global FIXEDPARAMETERS
% function that generates any of 20 peak types specified by number. 'shape'
% specifies the shape type of each peak in the signal: "peakshape" = 1-20.
% 1=Gaussian 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 9=exponential pulse, 10=up sigmoid,
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=Breit-Wigner-Fano (BWF) ,
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile; 21=triangular; 23=down sigmoid; 25=lognormal. "extra" is required
% for variable-shape peaks only.
switch shape,
    case 1
        p=gaussian(x,pos,wid);
    case 2
        p=lorentzian(x,pos,wid);
    case 3
        p=logistic(x,pos,wid);
    case 4
        p=pearson(x,pos,wid,extra);
    case 5
        p=expgaussian(x,pos,wid,extra);
    case 6
        p=gaussian(x,pos,wid);
    case 7
        p=lorentzian(x,pos,wid);
    case 8
        p=expgaussian(x,pos,wid,extra)';
    case 9
        p=exppulse(x,pos,wid);
    case 10
        p=upsigmoid(x,pos,wid);
    case 11
        p=gaussian(x,pos,wid);
    case 12
        p=lorentzian(x,pos,wid);
    case 13
        p=GL(x,pos,wid,extra);
    case 14
        p=BiGaussian(x,pos,wid,extra);
    case 15
        p=BWF(x,pos,wid,extra);
    case 16
        p=gaussian(x,pos,wid);
    case 17
        p=lorentzian(x,pos,wid);
    case 18
        p=explorentzian(x,pos,wid,extra)';
    case 19
        p=alphafunction(x,pos,wid);
    case 20
        p=voigt(x,pos,wid,extra);
    case 21
        p=triangular(x,pos,wid);
    case 23
        p=downsigmoid(x,pos,wid);
    case 25
        p=lognormal(x,pos,wid);
    case 26
        p=linslope(x,pos,wid);
    case 27
        p=d1gauss(x,pos,wid);
    case 28
        p=polynomial(x,coeff);
    case 30
        p=voigt(x,pos,pos,wid,extra);
    case 31
        p=expgaussian(x,pos,wid,-extra);
    case 32
        p=pearson(x,pos,wid,extra);
    case 33
        p=GL(x,pos,wid,extra);
    case 34
        % width(extra)=abs(FIXEDPARAMETERS(extra));
        p=voigt(x,pos,wid,extra);
    case 35
        p=GL(x,pos,wid,extra);
    case 36
        p=expgaussian(x,pos,wid,-extra);
    case 37
        p=pearson(x,pos,wid,extra);
    case 38
        p=explorentzian(x,pos,wid,-extra);
    case 40
        p=sine(x,pos,wid);
    case 41
        p=rectangle(x,pos,wid);
    case 42
        p=ngaussian(x,pos,wid,extra);
    case 43
        p=Gompertz(x,pos,wid,extra);
    case 44
        p=OneMinusExp(x,pos,wid);
    case 45
        p=FourPL(x,pos,wid,extra);
    case 46
        p=quadslope(x,pos,wid);
    case 47
        p=blackbody(x,pos);
    case 48
        p=exppulse(x,pos,wid);
    otherwise
end % switch


function [coef,RSquared,StdDevs,BootResults]=plotit(xi,yi,polyorder,datastyle,fitstyle)
% [coef,RSquared,StdDevs]=plotit(x,y,n,datastyle,fitstyle)
% plotit accepts data in the form of a single vector "y", a pair of vectors
% "x" and "y", or a 2xn or nx2 matrix with x in first row or column. If the
% input argument "polyorder" is supplied, it fits the data to a polynomial
% of order "polyorder", plots the data in red dots and the fit as a blue
% line, and displays the fit coefficients and R-squared in the upper left
% corner of the graph. Polyorder=1 for straight line, =2 for quadratic
% (parabola) etc. 
%   The 4th and 5th input arguments are strings defining the line and symbol
% style and color, in standard Matlab convention. The strings, in single
% quotes, are made from one element from any or all the following 3
% columns:
%  
%            b     blue          .     point              -     solid
%            g     green         o     circle             :     dotted
%            r     red           x     x-mark             -.    dashdot 
%            c     cyan          +     plus               --    dashed   
%            m     magenta       *     star             (none)  no line
%            y     yellow        s     square
%            k     black         d     diamond
%            w     white         v     triangle (down)
%                                ^     triangle (up)
%                                <     triangle (left)
%                                >     triangle (right)
%                                p     pentagram
%                                h     hexagram
% For example, plotit(x,y,3,'or','-g') plots the data as red circles and
% the fit as a green solid line.
%   If the 4th output argument (BootResults) is supplied, plotit computes
% coefficient error estimates by the bootstrap method and returns the
% results in the matrix "BootResults" (of size 5 x polyorder+1). You can
% change the number of bootstrap samples in line 93. 
% Tom O'Haver, Version 6.1; Version 6, March 2015, added optional specification of
% plot style and color. Version 5 added standard deviations of coefficients
% as 3rd output argument (moved bootstrap results to 4th output argument).
% Examples:
% x=[1 2 3 4 5];y=[0 2 3 3 5];
% plotit(y); % plot y only vs index number, no fit.
% plotit(x,y); % plot x vs y only, data in separate vectors
% plotit([x;y]); % plot data only, data in matrix with x in first row
% plotit([x;y]'); % plot data only, data in matrix with x in first column
% plotit(y,2); % plot y vs index and fit to second order polynomial
% plotit(x,y,2); % plot x vs y and fit t  second order polynomial
% plotit(x,y,3,'or'); % as above, plotting data as red circles
% plotit(x,y,3,'or','-g'); % data as red circles, fit as green line
% plotit([x;y],3);  % plot and fit data in matrix
% [coef, RSquared]=plotit([x;y],2) % return fit coefficients and r-squared
% [coef,RSquared,STD]=plotit([x;y],2) Returns vector standard devations
% of coefficients as 'STD'. (STD./coef computes relative standard deviation)
% [coef, RSquared,,StdDevs,BootResults]=plotit(x,y,polyorder) % computes
% coefficient error estimates by the bootstrap method and returns the
% results in the matrix "BootResults" 
%
% You can use plotit.m to linearize and plot nonlinear relationships, such as
% y = a*exp(b*x) ; [coeff,R2]=plotit(x,log(y),1);a=exp(coeff(2));b=coeff(1);
% y = a*ln(b*x) ; [coeff,R2]=plotit(log(x),y,1);a=coeff(1);b=log(coeff(2));
% y=ax^b ; [coeff,R2]=plotit(log(x),log(y),1);a=exp(coeff(2));b=coeff(1);
%
% Related functions (http://tinyurl.com/cey8rwh): 
% plotfita.m animates the bootstrap process for instructional purposes. 
% logplotfit.m plots and fits log(x) vs log(y), for data that follows a
%    power law relationship or that covers a very wide numerical range. 
% trypoly(x,y) fits a series of polynomials to the data x,y, for
%   polynomial orders 1 through length(x)-1, returns the R2 values in a
%   vector. Shows that, for any data, R2 -> 1 as order -> length(x)-1.
% trydatatrans(x,y,order) tries 8 different simple data transformations
%   on the data x,y, fits the transformed data to a polynomial of order
%   'polyorder', returns the R2 values in a vector. Displays results in
%   3 x 3 array of small plots.

% Copyright (c) 2015, Thomas C. O'Haver

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
% 
% process input arguments
NumTrials=100;
if nargin==1, % Single argument, might be y vector or x,y matrix
    polyorder=0;
    datasize=size(xi);
    if isvector(xi), % xi is vector, not matrix
        y=xi;
        x=1:length(y); % Use this only to create an x vector if needed
        xlabelstring='index number';
    else % x is a matrix, split it into x and y vectors
        if datasize(1)<datasize(2),xi=xi';end
        x=xi(:,1);
        y=xi(:,2);
        polyorder=0;
        xlabelstring='x';
    end
end

if nargin==2, % Two arguments, might be x,y vectors or matrix and polyorder
    datasize=size(xi);
    if isscalar(yi),   %  matrix and polyorder
        if datasize(1)<datasize(2),xi=xi';end
        if isvector(xi),
            x=1:length(xi);
            xlabelstring='index number';
            y=xi';
            polyorder=yi;
        else
            x=xi(:,1);
            y=xi(:,2);
            polyorder=yi;
            xlabelstring='x';
        end
    else %  x,y vectors
        x=xi;
        y=yi;
        polyorder=0;
        xlabelstring='x';
    end
end

if nargin>2, % Three arguments, must be x,y vectors and polyorder
    x=xi;
    y=yi;
    xlabelstring='x';
end

x=reshape(x,[1,length(x)]);
y=reshape(y,[1,length(y)]);

if length(x)==length(y),
    if polyorder>0, % Plot and fit the data
        % Compute the fit
        [coef,S]=polyfit(x,y,polyorder);
        StdDevs=sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';
        % Plot the data and the polynomial fit
        xx=linspace(min(x),max(x)); % 100 values of x over range of x values
        yhat=polyval(coef,xx); % Predicted values of y for each xx
        plot(x,y,'.r',xx,yhat,'-b')
        if nargin==4,
            plot(x,y,datastyle,xx,yhat,'-b')
        end
        if nargin==5,
            plot(x,y,datastyle,xx,yhat,fitstyle)
        end
        axis([min(x) max(x) min(y) max(y)]);
        xlabel(xlabelstring);ylabel('y')
        title(['plotit.m     Number of data points= ' num2str(length(x)) ] )
        % Compute the correlation coefficient and R-Squared
        if IsOctave,
            cc=corr(polyval(coef,x'),y');
            RSquared=cc.^2;
        else
            cc=corrcoef(polyval(coef,x),y);
            RSquared=cc(2).^2;
        end %   if IsOctave,
        % Label the graph with the fit information
        text(min(x),max(y)-.05.*(max(y)-min(y)),['   Polynomial Order of fit = ' num2str(polyorder)] );
        if polyorder==1;
            text(min(x),max(y)-.1.*(max(y)-min(y)),['   Slope = ' num2str(coef(1)) '   Intercept = ' num2str(coef(2))   ] );
        else
            text(min(x),max(y)-.1.*(max(y)-min(y)),['   Fit coefficients = ' num2str(coef)] );
        end
        text(min(x),max(y)-.15.*(max(y)-min(y)),['   Standard deviation = ' num2str(StdDevs)] );
        text(min(x),max(y)-.2.*(max(y)-min(y)),['   R-Squared = ' num2str(RSquared)] );
    else % Just plot the data, no fit.
        plot(x,y)
        axis([min(x) max(x) min(y) max(y)]);
        xlabel(xlabelstring);ylabel('y')
        title(['Number of data points= ' num2str(length(x)) ] )
        coef=[];
        RSquared=[];
    end
else
    disp('Plotit Error: x and y must be the same size.')
    sizex=size(x)
    sizey=size(y)
end
    
if nargout==4, % Compute the bootsrtap
    BootstrapResultsMatrix=zeros(polyorder+1,NumTrials);
    clear bx by
    for trial=1:NumTrials,
        n=1;
        bx=x;
        by=y;
        while n<length(x)-1,
            if rand>.5,
                bx(n)=x(n+1);
                by(n)=y(n+1);
            end % if rand>.5,
            n=n+1;
        end % while n<length(xx)-1,
        coef=polyfit(bx,by,polyorder);
        BootstrapResultsMatrix(:,trial)=coef;
    end % for trial=1:NumTrials,
    disp(' ')
    disp('Bootstrap Results');
    BootstrapMean=mean(real(BootstrapResultsMatrix(:,:)'));
    BootstrapSTD=std(BootstrapResultsMatrix(:,:)');
    BootstrapIQR=IQrange(BootstrapResultsMatrix(:,:)');
    PercentRSD=100.*BootstrapSTD./abs(BootstrapMean);
    PercentIQR=100.*BootstrapIQR./abs(BootstrapMean);
    BootstrapMean=BootstrapMean(1:polyorder+1);
    BootstrapSTD=BootstrapSTD(1:polyorder+1);
    BootstrapIQR=BootstrapIQR(1:polyorder+1);
    PercentRSD=PercentRSD(1:polyorder+1);
    PercentIQR=PercentIQR(1:polyorder+1);
        disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
        disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
        disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
        disp(['Percent RSD:    ', num2str(PercentRSD)])
        disp(['Percent IQR:    ', num2str(PercentIQR)])
    BootResults(1,:)=BootstrapMean;
    BootResults(2,:)= BootstrapSTD ;
    BootResults(3,:)=PercentRSD;
    BootResults(4,:)=BootstrapIQR;
    BootResults(5,:)=PercentIQR;
end % if NumArgOut==3,

function b=IQrange(a)
% b = IQrange(a) returns the interquartile range of the values in a.  For
%  vector input, b is the difference between the 75th and 25th percentiles
%  of a.  For matrix input, b is a row vector containing the interquartile
%  range of each column of a. (c) T. C, O'Haver, 2012
% Example:
%  a=randn(10000,5);
%  IQrange(a)
% Divide by 1.34896 to get an estimate of the standard deviation withput
% outliers
mina=min(a);
sizea=size(a);
NumCols=sizea(2);
b=zeros(size(a));
for n=1:NumCols,
    b(:,n)=a(:,n)-mina(n);
end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));

function isOctave = IsOctave()
% Returns true if this code is being executed by Octave.
% Returns false if this code is being executed by MATLAB, or any other MATLAB
% variant.
% 
%    usage: isOctave = IsOctave()
    persistent octaveVersionIsBuiltIn;
    if (isempty(octaveVersionIsBuiltIn))
        octaveVersionIsBuiltIn = (exist('OCTAVE_VERSION', 'builtin') == 5);
        % exist returns 5 to indicate a built-in function.
    end
    isOctave = octaveVersionIsBuiltIn;
    % If OCTAVE_VERSION is a built-in function, then we must be in Octave.
    % Since the result cannot change between function calls, it is cached in a
    % persistent variable.  isOctave cannot be a persistent variable, because it
    % is the return value of the function, so instead the persistent result must
    % be cached in a separate variable.
  
    
    function M=measurepeaks(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,plots)
% Peak detection and height and area measurement, similar syntax as
% findpeaksSG with the addition of a 7th input argument that controls
% plotting (0=no plots, 1=plots). The 3rd to 6th input arguments can be
% scalars or vectors with one entry for each segment.
% M=measurepeaks(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,plots)
% Locates and measures the positive peaks in a noisy x-y time series data.
% Detects peaks by looking for downward zero-crossings in the first
% derivative whose upward slopes exceed SlopeThreshold. Returns matrix (M)
% containing peak number, the absolute peak height, peak-valley difference,
% perpendicular drop area, and the tangent skim area of each peak.
% Arguments "SlopeThreshold", "AmpThreshold" and "SmoothWidth" control peak
% sensitivity of each segment. Higher values will neglect smaller features.
% "SmoothWidth" is a vector of the widths of the smooths applied before
% peak detection; larger values ignore narrow peaks. If smoothwidth=0, no
% smoothing is performed. "FitWidth" is a vector of the number of points
% around the top part of the peak that are taken for initial estimate of
% the peak center and width (minimum 3). If plots = 0, the function does
% not plot; if plots=1, the function plots the signal and numbers the peaks
% in figure 1 and plots each peak separately and marks the valley positions
% in an array of subplots in figure 2). Note: this function uses
% SmoothWidth only for peak detection; it performs measurements on the raw
% unsmoothed y data. If the data are noisy, it may be beneficial to smooth
% the y data yourself before calling measurepeaks.m, using any smooth
% function of your choice.
%
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver, 2016.  Version 1.1, December, 2016
%
% The script testmeasurepeaks.m will run all of the following examples with
% a 1-second pause between each:
%
% Example 1: sin(x).^2 has theoretical peaks at x=0.5pi, 1.5pi, 2.5pi,
% 3.5pi..., with peak heights of 1.0 and peak areas of pi/2 = 1.5708.
%  disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
%  x=[0:.01:20]';y=sin(x).^2;measurepeaks(x,y,0,0,5,5,1)
%
% Example 2: The built-in "humps" function has two peaks, no noise.
%  measurepeaks(0:.01:2,humps(0:.01:2),0,0,1,1,1)
%
% Example 3: Series of peaks that get progressively taller and wider.
%  x=[0:.01:5]';y=x.*sin(x.^2).^2;measurepeaks(x,y,0,0,5,5,1)
%
% Example 4: Like example 3, with random white noise added.
%  x=[0:.01:5]';y=.1.*randn(size(x))+x.*sin(-x.^2).^2;
%  measurepeaks(x,y,.001,.5,15,15,1)
%
% Example 5; Like example 3, with added rising baseline.
%  x=[0:.01:5]';y=x+x.*sin(x.^2).^2;measurepeaks(x,y,0,0,5,5,1)
%
% Example 6: Gaussian on linear baseline, theoretical area = 1.7725
%  x=[0:.1:10];y=2+x/10+exp(-(x-5).^2)+.01.*randn(size(x));
%  disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
%  disp(measurepeaks(x,y,0.0001,0.3,3,5,1))
%  Peak-valley and Tan skim are most accurate height and area measures.
%
% Example 7: Two overlapping Gaussians, zero baseline, theoretical area = 1.7725
%   x=[0:.05:10];y=exp(-(x-6).^2)+exp(-(x-3.5).^2)+.01.*randn(size(x));
%  disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
%  disp(measurepeaks(x,y,0.0001,0.5,8,6,1))
%  PeakMax and Perp drop are most accurate height and area measures.
% 
% Example 8: Narrow Gaussian peak on sloping linear baseline, followed by
% much broader peak, theoretical heights 1 and 1, areas 1.59 and 31.9.
% x=1:.2:150; 
% y=(150-x)./200+gaussian(x,30,1.5)+gaussian(x,90,30);
% measurepeaks(x,y,[0.001 .00005],[.6 .6],[5 120],[8 150],1)

if nargin==6, plots=1; disp('To supress plotting, set "plots" to 0.');end
if FitWidth<3,FitWidth=3;end
P=findpeaksSG(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
sizeP=size(P);
NumPeaks=sizeP(1);
clear PeakHeight
clear M
hold off
if plots,
    % Label the x-axis with the parameter values
    figure(1)
    clf
    title('measurepeaks.m function. Detected peaks are numbered. Peak table is printed in Command Window')
    xlabel('X');ylabel('Y')
    figure(2)
    clf
    hold off
end % if plots
% M=zeros(NumPeaks,6);
goodpeak=1;
for PeakNumber=1:NumPeaks,
    PeakPos=P(PeakNumber,2);
    PeakWidth=P(PeakNumber,4);
    PeakIndex=val2ind(x,PeakPos);
    x1=val2ind(x,PeakPos-1.3.*PeakWidth); % Left-most end of peak
    x2=val2ind(x,PeakPos+1.3.*PeakWidth);  % Right-most end of peak
    Halfway=val2ind(x,PeakPos); % x-value of center of peak
    miny=min(y(x1:x2)); % Minimum y value of that region
    yy(x1:x2)=y(x1:x2)-miny;  % Subtract minimum from that region
    % Calculate height, peak-valley, and perpendicular drop area
    PeakMax=y(PeakIndex);
    LeftSideMinimum=min(yy(x1:Halfway));
    RightSideMinimum=min(yy(Halfway:x2));
    LeftValleyIndex=x1+val2ind(yy(x1:Halfway),LeftSideMinimum);
    RightValleyIndex=PeakIndex+val2ind(yy(Halfway:x2),RightSideMinimum);
    if RightValleyIndex>length(x),RightValleyIndex=length(x);end
    PDarea=sum(y(LeftValleyIndex:RightValleyIndex)).*(x(2)-x(1)); % Perpendicular Drop area
    PVHeight=PeakMax-(interp1([x(LeftValleyIndex) x(RightValleyIndex)],[y(LeftValleyIndex) y(RightValleyIndex)],PeakPos));
    interpy=interp1([x(LeftValleyIndex) x(RightValleyIndex)],[y(LeftValleyIndex) y(RightValleyIndex)],x(LeftValleyIndex:RightValleyIndex));
    TSarea=sum((y(LeftValleyIndex:RightValleyIndex)-interpy).*(x(x1+1)-x(x1)));
    % [PeakNumber PeakPos PeakMax PVHeight]
    if PVHeight>AmpThreshold,
        M(goodpeak,:)=[goodpeak PeakPos  PeakMax PVHeight PDarea TSarea];
        if plots,
            % Plot individual peaks in a square subplot grid
            side=round(sqrt(NumPeaks)+.4); % Number of subplots on a side of square array
            if PeakNumber>side^2,
            else
            subplot(side,side,PeakNumber) % Set up side^2 square set of subplots
            end
            plot(x(x1:x2),y(x1:x2)) % Plot the individual peak
            title(['Peak ' num2str(goodpeak) ]) % Title the individual peak
            xlabel('X');ylabel('Y')
            % Mark valleys with magenta lines
            hold on
            if LeftValleyIndex<1,LeftValleyIndex=1;end
            plot([x(LeftValleyIndex) x(LeftValleyIndex)],[miny PeakMax],'m.-')
            if RightValleyIndex>length(x),RightValleyIndex=length(x);end
            plot([x(RightValleyIndex) x(RightValleyIndex)],[miny PeakMax],'m.-')
            plot([x(LeftValleyIndex) x(RightValleyIndex)],[interpy(1) interpy(length(interpy))],'c.-')
            plot(x(PeakIndex),PeakMax,'ro')
            hold off
        end % if plots
        goodpeak=goodpeak+1;
    end
end % for PeakNumber
figure(1)
if goodpeak-1>0,plot(x,y);text(M(:,2),M(:,3),num2str(M(:,1))),end

function P=findpeaksSG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Segmented peak finder, same syntax as findpeaksG except the 3rd to 6th
% input arguments can be vectors with one entry for each segment.
% function P=findpeaksSG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Locates and measures the positive peaks in a noisy x-y time series data.
% Detects peaks by looking for downward zero-crossings in the first
% derivative whose upward slopes exceed SlopeThreshold. Returns list (P)
% containing peak number and position, height, width, and area of each
% peak. Arguments "slopeThreshold", "ampThreshold" and "smoothwidth"
% control peak sensitivity of each segment. Higher values will neglect
% smaller features. "Smoothwidth" is a vector of the widths of the smooths
% applied before peak detection; larger values ignore narrow peaks. If
% smoothwidth=0, no smoothing is performed. "Peakgroup" is a vector of the
% number points around the top part of the peak that are taken for
% measurement. If Peakgroup=0 the local maximum is taken as the peak height
% and position. The argument "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar) If
%   smoothtype=2, triangular (2 passes of sliding-average) If smoothtype=3,
%   pseudo-Gaussian (3 passes of sliding-average)
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver, 2016.  Version 1, November, 2016
%
% Example: Find, measure, and plot 2 noisy peaks with very different widths
%   x=1:.2:100;
%   y=gaussian(x,20,1.5)+gaussian(x,80,30)+.02.*randn(size(x));
%   plot(x,y,'c.')
%   P=findpeaksSG(x,y,[0.001 .0001],[.2 .2],[5 10],[10 50],3)
%   text(P(:,2),P(:,3),num2str(P(:,1)))
%   disp(' ')
%   disp('           peak #    Position      Height     Width      Area')
%   disp(P)
%
% Related functions:
% findpeaksG.m, findvalleys.m, findpeaksL.m, findpeaksb.m, findpeaksb3.m,
% findpeaksplot.m, peakstats.m, findpeakSNR.m, findpeaksGSS.m,
% findpeaksLSS.m, findpeaksfit.m, findsteps.m, findsquarepulse.m, idpeaks.m

% Copyright (c) 2016 Thomas C. O'Haver
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
if nargin~=7;smoothtype=1;end  % smoothtype=1 if not specified in argument
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end
if smoothwidth<1;smoothwidth=1;end
if isscalar(AmpThreshold),AmpThreshold=AmpThreshold.*ones(size(SlopeThreshold));end
if isscalar(smoothwidth),smoothwidth=smoothwidth.*ones(size(SlopeThreshold));end
if isscalar(peakgroup),peakgroup=peakgroup.*ones(size(SlopeThreshold));end
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=SegmentedSmooth(deriv(y),smoothwidth,smoothtype);
else
    d=deriv(y);
end

P=[0 0 0 0 0];
vectorlength=length(y);
NumSegs=length(SlopeThreshold);
peak=1;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth-1,
    Seg=fix(1+NumSegs./(vectorlength./j));
    n=round(peakgroup(Seg)/2+1);
    % [j Seg]
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold(Seg), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpThreshold(Seg),  % if height of peak is larger than AmpThreshold
                xx=zeros(size(peakgroup(Seg)));yy=zeros(size(peakgroup(Seg)));
                for k=1:peakgroup(Seg), % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    % groupindex=groupindex
                    xx(k)=x(groupindex);
                    yy(k)=y(groupindex);
                end
                if peakgroup(Seg)>2,
                    [Height,Position,Width]=gaussfit(xx,yy);
                    PeakX=real(Position);   % Compute peak position and height of fitted parabola
                    PeakY=real(Height);
                    MeasuredWidth=real(Width);
                    % if the peak is too narrow for least-squares technique to work
                    % well, just use the max value of y in the sub-group of points near peak.
                else
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end
                % Construct matrix P. One row for each peak detected,
                % containing the peak number, peak position (x-value) and
                % peak height (y-value). If peak measurement fails and
                % results in NaN, or if the measured peak height is less
                % than AmpThreshold, skip this peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold(Seg),
                    % Skip this peak
                else % Otherwise count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth];
                    peak=peak+1; % Move on to next peak
                end
            end
        end
    end
end
% ----------------------------------------------------------------------

function d=deriv(a)
% First derivative of vector using 2-point central difference.
%  T. C. O'Haver, 1988.
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
    d(j)=(a(j+1)-a(j-1)) ./ 2;
end

function SmoothedSignal=SegmentedSmooth(y,smoothwidths,type,ends)
%   SegmentedSmooth(y,w,type,ends) divides y into a number of equal-length
%   segments defined by the length of the vector 'smoothwidths', then
%   smooths each segment with smooth of type 'type' and width defined by
%   the elements of vector 'smoothwidths'.
%    Version 1, November 2016.
%   The argument "type" determines the smooth type:
%     If type=1, rectangular (sliding-average or boxcar)
%     If type=2, triangular (2 passes of sliding-average)
%     If type=3, pseudo-Gaussian (3 passes of sliding-average)
%     If type=4, pseudo-Gaussian (4 passes of same sliding-average)
%     If type=5, multiple-width (4 passes of different sliding-averages)
%   The argument "ends" controls how the "ends" of the signal
%   (the first w/2 points and the last w/2 points) are handled.
%     If ends=0, the ends are zero. (In this mode the elapsed
%       time is the fastest and is independent of the smooth width).
%     If ends=1, the ends are smoothed with progressively
%       smaller smooths the closer to the end. (In this mode the
%       elapsed time increases with increasing smooth widths).
%   SegmentedSmooth(Y,w,type) smooths with ends=0.
%   SegmentedSmooth(Y,w) smooths with type=1 and ends=0.
%
%   Examples: 3-segment smooth of random white noise, smooth widths of
%   2,20, and 200.
%     x=1:10000;y=randn(size(x));
%     plot(x,SegmentedSmooth(y,[2 20 200],3,0))
%
%     20-segment smooth, odd smooth widths from 1 to 41:
%     plot(x,SegmentedSmooth(y,[1:2:41],3,0))

%   Copyright (c) 2012, Thomas C. O'Haver

if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
ly=length(y);
NumSegments=length(smoothwidths);
SegLength=round(ly./NumSegments);
SmoothSegment=zeros(ly,NumSegments);
SmoothedSignal=zeros(1,ly);
for Segment=1:NumSegments,
    SmoothSegment(:,Segment)=fastsmooth(y,smoothwidths(Segment),type,ends);
    startindex=(1+(Segment-1)*SegLength);
    endindix=startindex+SegLength-1;
    if endindix>ly,endindix=ly;end
    indexrange=startindex:endindix;
    SmoothedSignal(indexrange)=SmoothSegment(indexrange,Segment);
end

function SmoothY=fastsmooth(Y,w,type,ends)
% fastsmooth(Y,w,type,ends) smooths vector Y with smooth
%  of width w. Version 3.0, October 2016.
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar)
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
%   If type=4, pseudo-Gaussian (4 passes of same sliding-average)
%   If type=5, multiple-width (4 passes of different sliding-average)
% The argument "ends" controls how the "ends" of the signal
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively
%     smaller smooths the closer to the end. (In this mode the
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Examples:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
%
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%
% x=1:100;
% y=randn(size(x));
% plot(x,y,x,fastsmooth(y,5,3,1),'r')
% xlabel('Blue: white noise.    Red: smoothed white noise.')
%
% Copyright (c) 2012, Thomas C. O'Haver
%

switch type
    case 1
        SmoothY=sa(Y,w,ends);
    case 2
        SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
        SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
    case 4
        SmoothY=sa(sa(sa(sa(Y,w,ends),w,ends),w,ends),w,ends);
    case 5
        SmoothY=sa(sa(sa(sa(Y,round(1.6*w),ends),round(1.4*w),ends),round(1.2*w),ends),w,ends);
end

function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
    s(k+halfw-1)=SumPoints;
    SumPoints=SumPoints-Y(k);
    SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
        SmoothY(k)=mean(Y(1:(2*k-1)));
        SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
end

function [Height, Position, Width]=gaussfit(x,y)
% Converts y-axis to a log scale, fits a parabola
% (quadratic) to the (x,ln(y)) data, then calculates
% the position, width, and height of the
% Gaussian from the three coefficients of the
% quadratic fit.  This is accurate only if the data have
% no baseline offset (that is, trends to zero far off the
% peak) and if there are no zeros or negative values in y.
%
% Example 1: Simplest Gaussian data set
% [Height, Position, Width]=gaussfit([1 2 3],[1 2 1])
%    returns Height = 2, Position = 2, Width = 2
%
% Example 2: best fit to synthetic noisy Gaussian
% x=50:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
% [Height,Position,Width]=gaussfit(x,y)
%   returns [Height,Position,Width] clustered around 100,100,100.
%
% Example 3: plots data set as points and best-fit Gaussian as line
% x=[1 2 3 4 5];y=[1 2 2.5 2 1];
% [Height,Position,Width]=gaussfit(x,y);
% plot(x,y,'o',linspace(0,8),Height.*gaussian(linspace(0,8),Position,Width))
% Copyright (c) 2012, Thomas C. O'Haver
% To prevent problems from taking the log of zero or negative values,
% make the lowest value of y equal to 1% of the maximum value.
maxy=max(y);
for p=1:length(y),
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
logyyy=log(abs(y));
[coef,S,MU]=polyfit(x,logyyy,2);
c1=coef(3);c2=coef(2);c3=coef(1);
% Compute peak position and height or fitted parabola
Position=-((MU(2).*c2/(2*c3))-MU(1));
Height=exp(c1-c3*(c2/(2*c3))^2);
Width=norm(MU(2).*2.35703/(sqrt(2)*sqrt(-1*c3)));

    