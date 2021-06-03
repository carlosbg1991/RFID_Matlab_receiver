function [Positions,Heights,Widths,shapes,FittingError]=fitshape3(x,y,start)
% Simplified peak fitting example using fminsearch function, with 3
% iterated variables. Fits a model shape defined in the last line of this
% function to the data in x,y, with the first-guess values defined by
% vector 'start' = [position1 width1 shape1 position2 width2 shape2
% ...etc]. The shape must has three iterated variables (e.g. position,
% width, and shape). Number of peaks in the model is 1/3 the length of
% 'start'. Change the last line to the expression for your shape. 
% Version 1, T. C. O'Haver, toh@umd.edu, 2016
%
% Example 1: voigt profile (Variables are: Positions,Widths,shapes)
% x=0:.1:10; 
% y=voigt(x,5,1,.1) + .01*randn(1,length(x));
% [Positions,Heights,Widths,shapes,FittingError]=fitshape3(x,y,[5 1 1])
%
% Example 2: 4-parameter logistic (miny, maxy, inflection point (IP), slope)
% x=0:.1:20;
% miny=0;slope=5;ip=10;d=0;maxy=10;
% y = maxy+(miny-maxy)./(1+(x./ip).^slope)+.1*randn(size(x));
% start=[1 1 1]; 
% [Miny,Maxy,slope,IP,FittingError]=fitshape3(x,y,start)
%
% Example 3: Weighted sum of two 4-parameter logistics.
% miny=0;slope=5;ip=10;d=0;maxy=10;
% y1= maxy+(miny-maxy)./(1+(x./ip).^slope);
% miny2=0;maxy2=20;ip2=30;slope2=20;
% y2=maxy2+(miny-maxy2)./(1+(x./ip2).^slope2);
% y=y1+y2+.2*randn(size(x));
% start=[0 10 5 0 10 30];
% [Miny,Maxy,slope,IP,FittingError]=fitshape3(x,y,start)
% 
% Example 4: Exponentially modified Gaussian (EMG)
% For positive peaks, change line 75 to PEAKHEIGHTS = abs(A\y');
% function g = shapefunction(t,mu,s,lambda)
%  G=s.*lambda.*sqrt(pi/2).*exp(0.5.*(s.*lambda).^2-lambda.*(t-mu)).*erfc((1/sqrt(2)).*(s.*lambda-((t-mu)./s)));
%  g=G./max(G);
%
format compact
global PEAKHEIGHTS
% Call the fminsearch function with specified fitting function
FitResults=fminsearch(@(lambda)(fitfunction(lambda,x,y)),start);
NumPeaks=round(length(start)./3);

% Construct model from FitResults
for m=1:NumPeaks,
  A(m,:)=shapefunction(x,FitResults(3*m-2),FitResults(3*m-1),FitResults(3*m));
end
model=PEAKHEIGHTS'*A;

% Compute percent fitting error 
FittingError=100*norm(y-model)./(sqrt(length(x))*max(y));

% Plot data as dots, fitted model as line
plot(x,y,'.',x,model,'r-')
hold on;
for m=1:NumPeaks,
    plot(x,A(m,:)*PEAKHEIGHTS(m),'g');
end
hold off
for m=1:NumPeaks,
  Positions(m)=FitResults(3*m-2);
  Widths(m)=FitResults(3*m-1);
  shapes(m)=FitResults(3*m);
end
Heights=PEAKHEIGHTS';

% -------------------------------------------------
function err = fitfunction(lambda,t,y)
%  Fitting function for a signal consisting of overlapping peaks.
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = shapefunction(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');

function g = shapefunction(x,miny,slope,ip)
%  shapefunction(x,pos,wid) = peak peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Replace the last line with the expression for your shape,
% similar to the gaussian definition originally there.
g = 1+(miny-1)./(1+(x./ip).^slope); % Change this to change shape