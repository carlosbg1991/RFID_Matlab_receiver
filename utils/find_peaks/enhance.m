function Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth)
% Resolution enhancement function by derivative method. the
% arguments factor1 and factor 2 are 2nd and 4th derivative weighting
% factors. Larger values of factor1 and factor2 will reduce the 
% peak width but will cause artifacts in the baseline near 
% the peak.  Adjust the factors for the the best compromise. 
% Use minimum smooth width needed to reduce excess noise. 
% Use iSignal to adjust these factors 
% 
% Example:
%  x=[0:.01:18];
%  y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2);
%  Enhancedsignal=enhance(y,1000,1000000,3)
%  plot(x,y,x,Enhancedsignal,'r')
%
% Related functions: SegmentedSharpen.m

d2=secderiv(signal);  % Computes second derivative
d4=secderiv(d2);   % Computes fourth derivative
Enhancedsignal = signal-factor1.*tsmooth(d2,SmoothWidth)+...
factor2.*tsmooth(tsmooth(tsmooth(d4,SmoothWidth),SmoothWidth),SmoothWidth);

function d=secderiv(a)
% Second derivative of vector using 3-point central difference.
% Example: secderiv([0 1 2 3 4 3 2 1 0]) yields [ 0 0 0 0 -2 0 0 0 0]
%  T. C. O'Haver, 2006.
n=length(a);
d=zeros(size(a));
for j = 2:n-1;
  d(j)=a(j+1) - 2.*a(j) + a(j-1);
end
d(1)=d(2);
d(n)=d(n-1);

function s=tsmooth(Y,w)
%  tsmooth(Y,w) smooths vector Y by a triangular function of halfwidth w
%  T. C. O'Haver, 1988.
v=ones(1,w);v=conv(v,v);
S=conv(Y,v);
startpoint=(length(v) + 1)/2;
endpoint=length(Y)+startpoint-1;

s=S(startpoint:endpoint) ./ sum(v);
