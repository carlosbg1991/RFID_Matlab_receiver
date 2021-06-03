function err = fitblackbody(lambda,wavelength,y,handle)
%  Fitting function for a blackbody spectrum.
%  T. C. O'Haver, May 2008
global emissivity
radiance = 1.19111E+16*wavelength.^(-5)./(exp(14380000./(wavelength*lambda))-1);
emissivity = radiance'\y';
z = radiance*emissivity;
err = norm(z-y);