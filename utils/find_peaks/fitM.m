function err = fitM(lam,yobsd,Spectra,InstFun,StrayLight)
% Fitting function for broadened absorption of any number of components
% yobsd =  observed transmission spectrum (column vector)
% Sprecta = reference spectra for each component, one component/column
% InstFunction = Instrument function or slit function. (column vector)
% StrayLight = fractional stray light (scalar or column vector)
% Example (4-point spectra; True absorbance is 1.0): 
% options = optimset('TolX',0.000001);% absorbance=fminsearch(@(lambda)(fitM(lambda,[0.56529 0.38696 0.56529 0.73496]',[0.2 1 0.2 0.058824]',[1 0.5 0.0625 0.5]',.01)),1)
% yobsd, Spectra, and InstFunction must have same number of rows (wavelengths)
%  T. C. O'Haver, August 2006
global z
global c
A = StrayLight + (10 .^ -(Spectra*lam'));
fy=fft(A);
fa=fft(InstFun);
fy1=fy.*fa;
z=real(ifft(fy1))./sum(InstFun);
c = z\yobsd;
q = z*c;
err = norm(q-yobsd);