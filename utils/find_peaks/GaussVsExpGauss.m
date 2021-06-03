% Curve fitting with exponentially modified Gaussian shapes. Comparison of
% shapes 31 and 39. Both have the same shape but are parameterized
% differently. Shape 31 reports the width as the FWHM of the original
% Gaussian and shape 39 reports the standard deviation (sigma) of that
% Gaussian. Shape 31 reports the exponential factor an the number of data
% points and shape 39 reports the reciprocal of time constant in time
% units. You must have peakfit.m version 8.4, gaussian.m, expgaussian.m,
% findpeaksG.m, and halfwidth.m in the Matlab/Octave path. Optional
% high-quality plotting from Matlab using the downloadable PlotPub package.
figure(1)
% User changeable constants:
xincrement=.1;
x=1:xincrement:101;
xposition=10;
FWHM=5; % full width at half maximum
timeconstant=20; % time constant in number of points

% Generate Gaussian (y) and exponentially broadened Gaussian (ey)
y=gaussian(x,xposition,FWHM); % Original Gaussian
ey=expgaussian(x,xposition,FWHM,-timeconstant)'; % Gaussuan convoluted with exponential function
lengthy=length(y)

% Plot them using inbuilt plot function
figure(1);clf
tic
plot(x,y,x,ey)
title('Comparison of Gaussian (blue) and Exponentially broadened Gaussian (green)')
xlabel(['Exponential time constant = ' num2str(timeconstant)  ' points  (lambda= ' num2str(1/(xincrement.*timeconstant)) ')' ])
toc

% compare width
w1=halfwidth(x,y);
w2=halfwidth(x,ey);

% Compare areas
areay=trapz(x,y);
areaey=trapz(x,ey);

% Compare peak positons and heights by Gaussian fitting
ypeaks=findpeaksG(x,y,0.001,0.02,7,7,3);
x1=ypeaks(2);
y1=ypeaks(3);
eypeaks=findpeaksG(x,ey,0.001,0.001,7,7,3);
x2=eypeaks(2);
y2=eypeaks(3);

% Label peaks
% text(x1,.96.*y1,[' Position=' num2str(x1)  '   Height=' num2str(y1)  '   Halfwidth=' num2str(w1) '   Area='  num2str(areay) ],'FontSize',12)
% text(x2,.98.*y2,['  Position=' num2str(x2)  '   Height=' num2str(y2)  '   Halfwidth=' num2str(w2) '   Area='  num2str(areaey) ],'FontSize',12)

% figure(2)
% FitResults1=peakfit([x;ey],0,0,1,31,timeconstant,5,0);
% subplot(211)
% title('Curve fitting with exponential convoluted Gaussian (shape 31)')
% 
% figure(3)
% FitResults2=peakfit([x;ey],0,0,1,39,1/timeconstant,5,0);
% subplot(211)
% title('Curve fitting with exponential modified Gaussian expression (shape 39)')
% disp(' ')
% disp('Method    Position     Height    Halfwidth       Area      Exponential factor')
% disp(['Shape 31    ' num2str(x1)  '          ' num2str(y1)  '            ' num2str(w1) '           '  num2str(areay)  '         '  num2str(FitResults1(6))])
% disp(['Shape 39    ' num2str(x2)  '    ' num2str(y2)  '      ' num2str(w2) '      '  num2str(areaey)  '          '  num2str(FitResults2(6))])

% % Plot them using downloadable PlotPub library
% Lines 35-39 create a high-quality plot of the Figure(1) graphic using the
% "Plot" function as an alternative to the Matlab "plot" function; You must
% download the PlotPub library from
% http://masumhabib.com/blog/plotpub-publication-quality-graph-v2-0-released/
% After unzipping, place the contents of the "lib" folder in your path.
% Press Ctrl-S to save; Ctrl-P to print.
tic
plt = Plot(x,y,x,ey); % create the figure directly using data 
plt.BoxDim=[9 7]; % Set size of graphic in inches (for publication)
plt.XLabel = (['Exponential time constant = ' num2str(timeconstant)  ' points  (lambda= ' num2str(1/(xincrement.*timeconstant)) ')' ]);
plt.YLabel = 'Absorbance'; % ylabel 
plt.Title = 'Comparison of Gaussian (blue) and Exponentially broadened Gaussian (red)'; % plot title
toc