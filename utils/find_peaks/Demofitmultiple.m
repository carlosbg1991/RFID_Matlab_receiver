% Demos iterative fit to sets computer-generated noisy peaks of different
% types, knowing only the shape types and variable shape parameters of each
% peak. Iterated parameters are shape, height, position, and width of both
% peaks. Requires the fitmultiple and peakfunction.m functions.
% T. C. O'Haver, February 2014
format compact
global PEAKHEIGHTS AUTOZERO BIPOLAR
clear A
BIPOLAR=0;
AUTOZERO=0;
% Create a set of  peaks and add random noise
k=1; % sets interval between x-values; smaller k = more points in signal
t=-100:k:500;
noise=2;
% 'shape' specifies the shape type of each peak in the signal: "peakshape" = 1-20.
% 1=Gaussian 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 6=equal-width Gaussians; 7=Equal-width Lorentzians;
% 8=exponentionally broadened equal-width Gaussian, 9=exponential pulse,
% 10=sigmoid, 11=Fixed-width Gaussian, 12=Fixed-width Lorentzian;
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=BiLorentzian,
% 16=Fixed-position Gaussians; 17=Fixed-position Lorentzians;
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile.
shape=[2 1 2 2]; 
% The following vectors define the height, position, width, and "m" of each
% peak in the signal
Height=[200 100 50 50];
Position=[0 150 300 450];
Width=[100 100 100 100];  
m=[0 1 0 0]; % Needed only for variable-shape peaks, otherwize = 0.
% Compute the signal
NumPeaks=length(Height);
for j = 1:NumPeaks,
    A(j,:) = Height(j).*peakfunction(shape(j),t,Position(j),Width(j),m(j));
end 
y=sum(A);
y = y + noise.*randn(size(t));

% Perform an iterative fit using the FMINSEARCH function
% **The 'shape' and 'm' vectors are assumed to be known**
% Compute first guess (start) vector
n=length(t);
start=[];
startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(t);
for marker=1:NumPeaks,
    markx=startpos(marker);
    start=[start markx n/ (3.*NumPeaks)];
end % for marker
options = optimset('TolX',.0001,'TolFun',.00001,'Display','off','MaxFunEvals',1000);
% "tic" and "toc" measure how long the iterative fit takes.
tic
parameter=fminsearch(@(lambda)(fitmultiple(lambda,t,y,shape,m)),start,options);
toc
% Compute a model and plot it (blue line) along with 
% the original data (red points)
for j = 1:NumPeaks,
    A(j,:) = PEAKHEIGHTS(j).*peakfunction(shape(j),t,parameter(2*j-1),parameter(2*j),m(j));
end 
model=sum(A);
plot(t,y,'r.',t,A',t,model,'k')
title('Demofitmultiple.m: fitting a signal with multiple peak types')