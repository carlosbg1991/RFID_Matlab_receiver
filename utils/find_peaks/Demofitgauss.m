% Demos iterative gaussian fit to computer-generated single noisy gaussian peak.
% Variable parameters are height, position, and width of peak.
% Uses gaussian and fitgauss2 functions.
% T. C. O'Haver, May 2006, edited to work in Octave, August 2013
format compact
global c

% Create a single gaussian peak and add random noise
k=1; % sets interval between x-values; smaller k = more points in signal
t=[-100:k:400];
noise=50;
Height=100;
Position=100; 
Width=100;  
y = Height.*gaussian(t,Position,Width) + noise.*randn(size(t));

% Perform an iterative fit using the FMINSEARCH function
start=[length(t)/2 length(t)/2];    % generic initial guesses for peak position and width
options = optimset('TolX',0.1);  % Determines how close the model must fit the data
% "tic" and "toc" measure how long the iterative fit takes.
tic
% parameter=fminsearch('fitgauss2',start,options,t,y);
parameter=fminsearch(@(lambda)(fitgauss2(lambda,t,y)),start);
toc
Measured=[c parameter(1) parameter(2)] % The peak height is returned in the global variable c.

% Compute a model and plot it (blue line) along with 
% the original data (red points)
model=c.*gaussian(t,parameter(1),parameter(2));
plot(t,y,'r.',t,model,'b')