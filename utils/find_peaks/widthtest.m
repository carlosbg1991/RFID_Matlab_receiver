% WidthTest demonstrates the fact that constraining some of the peak
% parameters of a fitting model to fixed values, if those values are
% accurately known, improves that accuracy of measurement of the other
% parameters, even though it increases the fitting error. Requires the GL.m
% and peakfit.m functions (version 7.6 or later).

% Create simulated x,y data with two overlapping peaks (each a 50%
% Gaussian/Lorentzian blend) and add 1% white noise.
x=-1:.1:11;
True=[1            4            1            3       3.5759           50
      2            6          0.5            3        1.788           50];
y=GL(x,4,3,50)+.5*GL(x,6,3,50)+.01*randn(size(x));
% Fit these x,y data with three different 2-peak models and plot results in
% Figures 1,2, and 3 and compute the relative percent parameter errors.
figure(1)
[FitResults13,FitError13]=peakfit([x;y],0,0,2,13,50,10);
subplot(2,1,1);title('Fixed shape factor (50% Gaussian) and variable widths')
Accuracy13=100.*mean(abs((True(:,2:5)-FitResults13(:,2:5))./True(:,2:5)));

figure(2);
[FitResults35,FitError35]=peakfit([x;y],0,0,2,35,50,10,0,0,[3 3]);
subplot(2,1,1);title('Fixed shape factor and fixed widths')
Accuracy35=100.*mean(abs((True(:,2:5)-FitResults35(:,2:5))./True(:,2:5)));

figure(3);
[FitResults33,FitError33]=peakfit([x;y],0,0,2,33,50,10,0,0,[3 3]);
subplot(2,1,1);title('Variable shape factor and widths')
Accuracy33=100.*mean(abs((True(:,2:6)-FitResults33(:,2:6))./True(:,2:6)));

% Display a table showing the percent fitting errors and the relative percent
% errors of peak position, height, width, and area for each model.
disp('Relative percent errors:')
disp('Fitting error  Position error  Height error  Width error  Area error  Shape error')
disp('Unconstrained shape factor and widths: shape 33, Figure 1')
disp([FitError33(1) Accuracy33])

disp('Fixed shape factor (50% Gaussian) and variable widths: Shape 13, Figure 2')
disp([FitError13(1) Accuracy13])

disp('Fixed shape factor and fixed widths [3 3]: Shape 35, Figure 3')
disp([FitError35(1) Accuracy35])

disp('Conclusion: the more constraints, the lower the parameter errors, IF the')
disp('constraints are accurate, but the fitting error will be slightly greater.')
disp('Run widthtest again to see the effect of a different noise sample.')