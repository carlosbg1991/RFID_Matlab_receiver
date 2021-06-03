function [a,b,FittingError,R2]=fitshape1(x,y,start)
% Simplified peak fitting example using fminsearch function. Fits a model
% shape defined the the shape function (line 60) to the data in x,y, with
% the first-guess values defined by vector 'start', in the form [position1
% position2, etc]. The shape must mave one iterated variable (e.g.
% position). Number of peaks in the model equals the length of 'start'.
% Change the last line to the expression for your shape. 
% Version 1, T. C. O'Haver, toh@umd.edu, 2016
%
% Example 1: Compound Interest
% NOTE: Change the last line to the expression for compound interest. 
% x=1:66; % time in years
% rate=.07; % Annual rate of return
% Noise=.5;
% y=start.*(1+rate).^x; %  Compound Interest formula
% yy=y+Noise.*y.*randn(size(x)); % Add proportional noise
% rate=fitshape1(x,yy,0)
%
% Example 2: single blackbody source
% NOTE: Change the last line to the expression for blackbody emission. 
% wavelength=[200 250 300 350 400 450 500 600]; % Wavelength in nm
% radiance = [.007 .011 .02 .04 .07 .08 .1 .1]; % Measured radiance in Watts/cm2/sr
% [Temperature,Emissivity,FittingError]=fitshape1(wavelength,radiance,3000)
% plot(x,yy,x,Heights*exp(Positions.*x))plot(x,yy,x,Heights*exp(Positions.*x))
%
% Example 3: Simulation of blackbody source with TWO temperature zones
% NOTE: Change the last line to the expression for blackbody emission. 
% w=200:1000;
% t1=5555;t2=3333;
% e1=.1;e2=1;
% y=blackbody(w,t1,e1)+blackbody(w,t2,e2);
% y=y+sqrt(y).*.01.*randn(size(w)); % Add photon noise
% [Temperatures,Emissivities,FittingError]=fitshape1(w,y,[3000 4000]);
%
format compact
global B
% Call the fminsearch function with specified fitting function
options = optimset('TolX',.000001,'TolFun',.00001,'Display','off','MaxFunEvals',1000);
FitResults=fminsearch(@(lambda)(fitfunction(lambda,x,y)),start,options);
NumPeaks=round(length(start));

% Construct model from FitResults
for m=1:NumPeaks,
  A(m,:)=shapefunction(x,FitResults(m));
end
model=B'*A;

% Compute percent fitting error 
FittingError=100*norm(y-model)./(sqrt(length(x))*max(y));
% Compute R2
SStot=sum((y-mean(y)).^2);
SSres=sum((y-model).^2);
R2=1-(SSres./SStot);

% Plot data as dots, fitted model as line
plot(x,y,'.',x,model,'r-')
hold on;
for m=1:NumPeaks,
    plot(x,A(m,:)*B(m),'g');
end
hold off
title('fitshape1')
for m=1:NumPeaks,
  a(m)=FitResults(m);
end
b=B';

% -----------------------------------------------------
function err = fitfunction(lambda,t,y)
%  Fitting function for a signal consisting of overlapping peaks.
global B
A = zeros(length(t),round(length(lambda)));
for j = 1:length(lambda),
    A(:,j) = shapefunction(t,lambda(j))';
end
B = A\y';
z = A*B;
err = norm(z-y');
% ----------------------------------------------------
function g = shapefunction(x,a)
%  shapefunction(x,a)
%  x may be scalar, vector, or matrix, a is scalar
%  T. C. O'Haver, 2016
% Some examples of shape functions to Copy and Paste into last line:
%
% Fixed-width GAUSSIAN peak, w=full width at half maximum (constant)
% g = exp(-((x-a)./(0.60056120439323.*w)) .^2);
%
% BLACKBODY emission expression: x=wavelength; a=temperature
% g = 1.19111E+16*x.^(-5)./(exp(14380000./(x*a))-1);
%
% EXPONENTIAL GROWTH (a > 0) OR DECAY (a < 0) EXPRESSION:
% g = exp(a.*x);
%
% COMPOUND INTEREST (a=interest rate)
% g=(1+a).^x;
%
g=(1+a).^x; % <<< CHANGE THIS LINE to the expression you want to fit
