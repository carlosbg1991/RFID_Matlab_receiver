function [Positions,Heights,Widths,FittingError]=fitshape2(x,y,start)
% Simplified peak fitting example using fminsearch function. Fits a model
% shape defined the the shape function (last line of this function) to the
% data in x,y, with the first-guess values defined by vector 'start' =
% [position1 width1 position2 width2 ...etc]. The shape must have three
% iterated variables (.eg. position, width, shape). Number of peaks in the
% model is 1/3 the length of 'start'. Change the last line to the
% expression for your shape. 
% Version 1, T. C. O'Haver, toh@umd.edu, 2016
%
% Example of use with a 2-peak Gaussian signal with noise:
% x=0:.1:10; 
% y=exp(-(x-5).^2)+2*exp(-(x-3).^2) + .1*randn(1,length(x));
% [Positions,Heights,Widths,FittingError]=fitshape(x,y,[2 1 6 1])
%
% Example of use with a 3-peak Gaussian signal with noise:
% x=0:.05:10;
% y=exp(-(x-5).^2)+2*exp(-(x-3).^2)+exp(-(x-7).^2) + .05*randn(1,length(x));
% [Positions,Heights,Widths,FittingError]=fitshape(x,y,[2 1 5 1 8 1])

format compact
global PEAKHEIGHTS
% Call the fminsearch function with specified fitting function
FitResults=fminsearch(@(lambda)(fitfunction(lambda,x,y)),start);
NumPeaks=round(length(start)./2);

% Construct model from FitResults
for m=1:NumPeaks,
  A(m,:)=shapefunction(x,FitResults(2*m-1),FitResults(2*m));
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
  Positions(m)=FitResults(2*m-1);
  Widths(m)=FitResults(2*m);
end
Heights=PEAKHEIGHTS';

% -------------------------------------------------
function err = fitfunction(lambda,t,y)
%  Fitting function for a signal consisting of overlapping peaks.
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = shapefunction(t,lambda(2*j-1),lambda(2*j))';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');

function g = shapefunction(x,pos,wid)
%  shapefunction(x,pos,wid) = peak peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Replace the last line with the expression for your shape,
% similar to the gaussian definition originally there.
g = exp(-((x-pos)./(0.6005612.*wid)) .^2); % Change this to change shape