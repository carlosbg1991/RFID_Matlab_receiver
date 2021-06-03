x=[0:.1:25];y=alphafunction(x,2,8)+.05*randn(size(x));
disp('Fitting errors')

% 'Two Gaussians'
figure(1)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,2,1,0,10,0,0,0);
disp(['Two Gaussians:          ' num2str(MeanFitError) ])

% 'Three Gaussians'
figure(2)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,3,1,0,10,[9 2 11 3 14 6],0,0); 
disp(['Three Gaussians:        ' num2str(MeanFitError) ])

% 'Variable Exponentially modified Gaussian'
figure(3)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,31,35,10,0,0,0);
disp([''Variable Exp. modified Gaussian: ' num2str(MeanFitError) ])

% 'BiGaussian'
figure(4)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,14,13,10,0,0,0);
disp(['BiGaussian:             ' num2str(MeanFitError) ])

% 'Exponential pulse'
figure(5)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,9,0,20,0,0,0);
disp(['Exponential pulse:      ' num2str(MeanFitError) ])

% 'Alpha function'
figure(6)
[FitResults,MeanFitError]=peakfit([x' y'],0,0,1,19,0,20,0,0,0);
disp(['Alpha function:         ' num2str(MeanFitError) ])

% 'Lognormal'
figure(7)
 [FitResults,MeanFitError]=peakfit([x' y'],0,0,1,25,35,10,0,0,0);
disp([''Lognormal':         ' num2str(MeanFitError) ])

