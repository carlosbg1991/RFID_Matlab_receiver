% Demonstration of exponentially modified Gaussian peaks
% Requires the following functions in the path: gaussian.m, expgaussian.m,
% halfwidth, val2ind.m, and peakfit.m
format compact; format short g
clear
increment=.05; % difference bewtween adjacent x-axis values
x=1:increment:50; % independent variable axis (e.g. time)
Hg=1; % peak height of original Gaussian
pg=15; % position of original Gaussian
wg=3; % width (FWHM) of original Gaussian
ExpFactor=40.1; % Exponential factor, in number of data points
% Create Gaussian, height Hg, centered at pg, FWHM=wg
g=Hg.*gaussian(x,pg,wg); 
% expgaussian, generated from gaussian of height Hg, centered at pg, FWHM=wg
ge=Hg.*expgaussian(x,pg,wg,-ExpFactor)'; 
figure(1)
clf
plot(x,g,x,ge)
axis([5 35 0 1])
title('Blue = Gaussian peak "g"    Red = Result of exponential broadening of "g" ')
% Determine peak parameters of both peaks
% Both peak areas are equal
disp('"g" = Gaussian  "ge" = exponentially braodened Gaussian ')
AreaOfg=sum(g).*increment
AreaOfge=sum(ge).*increment
% Peak position of ge
maxge=max(ge); % Max height of ge
val2ind(ge,maxge); % index (point  number) where ge equals maxge
PeakPositionOfge=x(val2ind(ge,maxge)) % Peak position of ge, in x units, is slightly higher than g
HeightOFge=maxge; 
% Peak widths (FWHM)
WidthOfg=halfwidth(x,g) 
WidthOfge=halfwidth(x,ge) 
% Plot of both peaks normalized and peak maxima coinciding
figure(2) 
clf
shift=PeakPositionOfge-pg % Shift in peak position
plot(x,g,x-shift,ge./maxge)
axis([5 35 0 1])
xlabel('x')
ylabel('y')
title('Plot of both peaks normalized and peak maxima coinciding')

% Iterative curve fitting of ge
figure(3)
clf
PeakShape=31;
[FitResults31,GOF31]=peakfit([x;ge],20,30,1,PeakShape,0,10);
subplot(2,1,1)
title('Curve fitting with exponentially broadened Gaussian (shape 31)')
disp('Curve fitting with exponentially broadened Gaussian (shape 31)')
disp(table(FitResults31(2),FitResults31(3),FitResults31(4),FitResults31(5),FitResults31(6),'VariableNames',{'Position' 'Height' 'FWHM' 'Area' 'ExpFactor'}))
disp(' ')
figure(4)
clf
PeakShape=39;
[FitResults39,GOF39]=peakfit([x;ge],20,30,1,PeakShape,0,10);
subplot(2,1,1)
title('Curve fitting with exponentially modified Gaussian (EMG) function (shape 39)')
disp('Curve fitting with exponentially modified Gaussian (EMG) function (shape 39)')
disp(table(FitResults39(2),FitResults39(3),FitResults39(4),FitResults39(5),FitResults39(6),'VariableNames',{'Position' 'Height' 'Sigma' 'Area' 'Lambda'}))
disp('Note that Sigma=FWHM/2.355 and that Lambda=ExpFactor.*increment.')
disp(' ')
disp('Conclusions:')
disp('Exponential broadening of a Gaussian peak shifts the peak to larger')
disp('x values, decreases the peak height, and increases the peak width,')
disp('but the area of the peak remains the same. Such broadened peaks can')
disp('be measured by the curve fitting function peakfit.m in two ways:')
disp('using shape number 31 will re-create the parameters of the original')
disp('Gaussian before the broadening, whereas using shape number 39 fits')
disp('the peak with the expression EMG = s.*lambda.*sqrt(pi/2).*exp(0.5.*')
disp('(s.*lambda).^2-lambda.*(t-mu)).*erfc((1/sqrt(2)).*(s.*lambda-((t-mu)./s)))')
