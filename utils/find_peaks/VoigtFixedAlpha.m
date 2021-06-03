function VoigtFixedAlpha
global PEAKHEIGHTS AUTOZERO
format short g
format compact
warning off all
% Signal definition (change as desired)
x=0:.01:10;
position1=3; % Peak position of first peak
DopWid1=.2; % Doppler width of first peak
alpha1=1.7; % alpha of first peak
position2=7;% Peak position of second peak
DopWid2=.3; % Doppler width of second peak
alpha2=1.7; % alpha of second peak
noise=0.01; % Standard deviation of random noise added to signal
y=voigt(x,position1,DopWid1,alpha1)+.5.*voigt(x,position2,DopWid2,alpha2)+noise.*randn(size(x));

% Curve fitting
NumPeaks=2;
AUTOZERO=0;
n=length(x);
newstart=[3 0.5 7 .9];
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+randn/100);
        newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
    end

options = optimset('TolX',.0001,'Display','off','MaxFunEvals',200 );
tic
TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,x,y,alpha1)),newstart);
toc
disp('       positon1      DopWid1     position2     DopWid2 ')
disp(TrialParameters)
alphas=[alpha1 alpha2];
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
    A(m,:)=voigt(x,TrialParameters(2*m-1),TrialParameters(2*m),alpha1);        
end
% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
end
alphas=alphas
% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(y-model)./(sqrt(n)*max(y))
  
  plot(x,y,'.k',x,model,'r')
  title('Voigt Test, fixed alpha')
  
  % ----------------------------------------------------------------------
function start=calcstart(xx,NumPeaks,xoffset)
  n=max(xx)-min(xx);
  start=[];
  startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
  for marker=1:NumPeaks,
      markx=startpos(marker)+ xoffset;
      start=[start markx n/ (3.*NumPeaks)];
  end % for marker
% ----------------------------------------------------------------------
function err = fitvoigt(lambda,t,y,shapeconstant)
% Fitting functions for Voigt profile function
% T. C. O'Haver (toh@umd.edu), 2013.
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g=voigt(xx,pos,gD,alpha)
% Voigt profile function. xx is the independent variable (energy,
% wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
% shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
% Based on Chong Tao's "Voigt lineshape spectrum simulation", 
% File ID: #26707
% alpha=alpha
gL=alpha.*gD;
gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2);
x = gL/gV;
y = abs(xx-pos)/gV;
g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
g=g./max(g);