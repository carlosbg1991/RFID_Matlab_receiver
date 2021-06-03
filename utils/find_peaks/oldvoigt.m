function g=oldvoigt(xx,pos,gD,alpha)
% Previous version of voigt.m used in peakfit.m 9.1 and earlier.
% Voigt profile function. xx is the independent variable (energy,
% wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
% shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
% Based on Chong Tao's "Voigt lineshape spectrum simulation", File ID:
% #26707
% alpha=alpha
gD=gD/2;
gL=alpha.*gD;
gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2);
ratio = gL/gV;
y = abs(xx-pos)./gV;
g = 1/(2*gV*(1.065 + 0.447*ratio + 0.058*ratio^2))*((1-ratio)*exp(-0.693.*y.^2) + (ratio./(1+y.^2)) + 0.016*(1-ratio)*ratio*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
g=g./max(g);