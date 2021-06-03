function results=cls2(xx,yy,NumPeaks,peakshape,Positions,Widths,extra)
% The function cls.m computes such a model ('modelmatrix') consisting of
% the sum of any number of peaks of known shape, width, and position, but
% of unknown height, and fits the model to the x,y data set. The syntax is 
% 
% results=cls2(x,y,NumPeaks,PeakShape,Positions,Widths,extra)
% 
% where x and y are the vectors of measured data, 'NumPeaks' is the number
% of peaks, 'PeakShape' is the peak shape number (1=Gaussian, 2=Lorentzian,
% 3=logistic, 4=Pearson, 5=exponentionally broadened Gaussian;
% 6=equal-width Gaussians; 7=Equal-width Lorentzians; 8=exponentionally
% broadened equal-width Gaussian, 9=exponential pulse, 10=sigmoid,
% 11=Fixed-width Gaussian, 12=Fixed-width Lorentzian;
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=BiLorentzian),
% 'Positions' is the vector of peak positions on the x axis (one entry per
% peak), 'Widths' is the vector of peak widths in x units (one entry per
% peak), and 'extra' is the additional shape parameter required by the
% exponentionally broadened, Pearson, Gaussain/Lorentzian blend, BiGaussain
% and BiLorentzian shapes. Cls2.m returns a vector containing the
% background B and measured peak heights h for each peak , e.g. [B h1 h2
% h3...]. Version 2, February 2016 toh@umd.edu
%
background=ones(size(yy));
results=yy/[background;modelmatrix(xx,NumPeaks,peakshape,Positions,Widths,extra)];

function A=modelmatrix(xx,NumPeaks,peakshape,Positions,Widths,extra)
% Specifies the peak shape of the model: "peakshape" = 1-15. (1=Gaussian
% (default), 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 6=equal-width Gaussians; 7=Equal-width Lorentzians;
% 8=exponentionally broadened equal-width Gaussian, 9=exponential pulse,
% 10=sigmoid, 11=Fixed-width Gaussian, 12=Fixed-width Lorentzian;
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=BiLorentzian
% Extra is needed only for shapes 4, 5, 8, 13, 14, and 15.
% FIXEDWIDTH is needed only for shapes 11 ands 12
%
A=zeros(NumPeaks,length(xx));
for m=1:NumPeaks,
    switch peakshape
        case 1
            A(m,:)=gaussian(xx,Positions(m),Widths(m));
        case 2
            A(m,:)=lorentzian(xx,Positions(m),Widths(m));
        case 3
            A(m,:)=logistic(xx,Positions(m),Widths(m));
        case 4
            A(m,:)=pearson(xx,Positions(m),Widths(m),extra);
        case 5
            A(m,:)=expgaussian(xx,Positions(m),Widths(m),-extra)';
        case 6
            A(m,:)=gaussian(xx,Positions(m),Widths(m));
        case 7
            A(m,:)=lorentzian(xx,Positions(m),Widths(m));
        case 8
            A(m,:)=expgaussian(xx,Positions(m),Widths(m),-extra)';
        case 9
            A(m,:)=exppulse(xx,Positions(m),Widths(m));
        case 10
            A(m,:)=sigmoid(xx,Positions(m),Widths(m));
        case 11
            A(m,:)=gaussian(xx,Positions(m),FIXEDWIDTH);
        case 12
            A(m,:)=lorentzian(xx,Positions(m),FIXEDWIDTH);
        case 13
            A(m,:)=GL(xx,Positions(m),Widths(m),extra);
        case 14
            A(m,:)=BiGaussian(xx,Positions(m),Widths(m),extra);
        case 15
            A(m,:)=BiLorentzian(xx,Positions(m),Widths(m),extra);
    end % switch
end % for
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);
% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function.
% g = pearson7(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
%  T. C. O'Haver, 1990
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2006
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) convolutes y by an exponential decay of time constant t
% by multiplying Fourier transforms and inverse transforming the result.
ly=length(y);
ey=[zeros(size(y));y;zeros(size(y))];
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(ly+2:length(ybz)-ly+1);
% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form
% Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* [p>0];
g = p';
% ----------------------------------------------------------------------
function g=sigmoid(x,t1,t2)
g=1./(1 + exp((t1 - x)./t2))';
% ----------------------------------------------------------------------
function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend. m = percent Gaussian character
% pos=position; wid=half-width
% m = percent Gaussian character.
%  T. C. O'Haver, 2012
g=((m/100)*gaussian(x,pos,wid)+(1-(m/100))*lorentzian(x,pos,wid))/2;
% ----------------------------------------------------------------------
function g = BiGaussian(x,pos,wid,m)
% BiGaussian (different widths on leading edge and trailing edge).
% pos=position; wid=width
% m determines shape; symmetrical if m=50.
%  T. C. O'Haver, 2012
lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=gaussian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,wid*(1-m/100));
% ----------------------------------------------------------------------
function g = BiLorentzian(x,pos,wid,m)
% BiLorentzian (different widths on leading edge and trailing edge).
% pos=position; wid=width
% m determines shape; symmetrical if m=50.
%  T. C. O'Haver, 2012
lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=lorentzian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=lorentzian(x(hx+1:lx),pos,wid*(1-m/100));
