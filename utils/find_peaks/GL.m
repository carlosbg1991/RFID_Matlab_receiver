function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend.
% pos=position; wid=half-width;
% m = percent Gaussian character.
%  T. C. O'Haver, 2012
g=2*((m/100)*gaussian(x,pos,wid)+(1-(m/100))*lorentzian(x,pos,wid))/2;