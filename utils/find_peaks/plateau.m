function g = plateau(x,pos,wid,sharp)
%  plateau(x,pos,wid,sharp) = flattened plateau centered on x=pos,
%  half-width=wid, sharp=risetime (sharpness)
%  x may be scalar or vector; pos, wid, and n are scalar
%  T. C. O'Haver, 2016
% Example
% x=[-5:.1:15];y=plateau(x,5,3,.4);plot(x,y);
% Independent check on half-width: halfwidth(x,y)

% up step sigmoid
g1=1/2 + 1/2* erf(real((x-pos+wid/2)/sqrt(2*sharp))); 
 % down step sigmoid
g2=.5-.5*erf(real((x-pos-wid/2)/sqrt(2*sharp)));
g=g1.*g2;