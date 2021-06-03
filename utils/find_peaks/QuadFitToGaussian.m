% Method of one-step least-squares fit to a single gaussian.
% Generates a noisy Gaussian function, converts the
%  y-axis to a log scale, fits a parabola
% (quadratic) to the (x,log(y)) data, then calculates
% the position, width, and height of the original
% Gaussian from the three coefficients of the
% quadratic fit.  This works only if the Gaussian has
% no baseline offset (that is, goes to zero far off the
% peak) and if there are no negative values in y.
% Tom O'Haver, toh@umd.edu, 2008

format compact
% Simulation section
% Generates a noisy Gaussian function
x=[25:31:200];
Height=100;
Width=100;
Position=100;
Noise=1;
% Newly-created example with random noise added
y=Height.*gaussian(x,Position,Width)+Noise.*randn(size(x));
% Example in ResultiingDataSet.txt on CurveFitting.html
y=[21.56 60.3 93.16 92.27 51.71 15.650];
% Analysis section;
z=log(abs(y));
coef=polyfit(x,z,2);
x2=[-50:300];
fit=polyval(coef,x2);
a=coef(3);
b=coef(2);
c=coef(1);
MeasuredHeight=exp(a-c*(b/(2*c))^2);
MeasuredPosition=-b/(2*c);
MeasuredWidth=2.35703/(sqrt(2)*sqrt(-c));
HeightResults=[Height MeasuredHeight];
PositionResults=[Position MeasuredPosition];
WidthResults=[Width MeasuredWidth];

% Plotting section
subplot(2,2,1) % Upper left segment
plot(x,y,'ro')
axis([-50,300,0,120]);
ylabel('y')
title('Original data')

subplot(2,2,2) % Upper right segment
plot(x,z,'ro')
axis([-50,300,2,5]);
ylabel('ln(y)')
title('Ln of original data')

subplot(2,2,3) % Lower left segment
plot(x,z,'ro',x2,fit)
axis([-50,300,2,5]);
ylabel('ln(y)')
title('Ln of original data with quadratic fit')
xlabel(['a = '  num2str(coef(3))  '    b = '   num2str(coef(2))  '    c = '  num2str(coef(1))] )

subplot(2,2,4) % Lower right segment
plot(x,y,'ro',x2,Height.*gaussian(x2,MeasuredPosition,MeasuredWidth),'b')
axis([-50,300,0,120]);
xlabel(['Height = '  num2str(MeasuredHeight)  '    Position = '   num2str(MeasuredPosition)  '    Width = '  num2str(MeasuredWidth)] )
ylabel('y')
title('Original data with computed Gaussian fit')
disp('  Height Error   Position Error  Width Error')
PercentErrors=[100*(MeasuredHeight-Height)/Height 100*(MeasuredPosition-Position)/Position 100*(MeasuredWidth-Width)/Width];
disp(PercentErrors)