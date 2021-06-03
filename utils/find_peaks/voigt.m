function v=voigt(x,pos,Gausswidth,voigtalpha)
% Unit height Voigt profile function. x is the independent variable
% (energy, wavelength, etc), Gausswidth is the Gaussian(Doppler) width,
% and voigtalpha is ratio of the Gausswidth to the LorentzWidth (pressure
% width). Version 3, August, 2019
LorentzWidth=Gausswidth.*voigtalpha;
if LorentzWidth<0, LorentzWidth=-LorentzWidth;end
if Gausswidth<0, Gausswidth=-Gausswidth;end
dx=x(2)-x(1); % x increment
ex=[x-max(x)-dx x x+max(x)]; % Extended x
gau=gaussian(ex,0,Gausswidth);
lor=lorentzian(ex,pos,LorentzWidth);
VoigtConv=ifft(fft(gau).*fft(lor))./sum(lor);
g=VoigtConv./max(VoigtConv);
oex=ex-max(x);
outrange=val2ind(oex,0):val2ind(oex,max(x))-dx;
v=g(outrange+1);