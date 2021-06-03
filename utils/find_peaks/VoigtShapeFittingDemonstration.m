format compact
format short g
x=0:.1:12;
pos=6;
G=1; % Gaussian width
L=1; % Lorentzian width
alpha=L./G % calculated alpha
% Create a single Voigt peak with known Gaussian and Lorentzian widths
v=voigt(x,pos,G,alpha);
clf
%  plot(x,v,'o')
disp('If the alpha is known amd the same from peak to peak, use shape 20')
NumPeaks=1;
PeakShape=20; % Voigt with fixed alpha
disp('          peak#       pos         height       width        area')
peakfit([x;v],0,0,NumPeaks,PeakShape,alpha)
disp(' ')
disp('The full width at half maximum (FWHM) of the Voigt peak can be approximated')
disp('from the widths of the known Gaussian and Lorentzian components')
disp('CalculatedFWHM = (0.5346*L + sqrt(0.2166*L.^2 + G.^2))')
CalculatedFWHM = (0.5346*L + sqrt(0.2166*L.^2 + G.^2))
MeasuredFWHM=halfwidth(x,v)
L=G.*alpha
area=trapz(x,v)
w=MeasuredFWHM;
disp('The widths of the Gaussian and Lorentzian components and be calulated knowing')
disp('one of those widths and the measured half maximum (FWHM) of the Voigt peak')
disp('CalculatedLorWidth=(5000*(2673*w - sqrt(3)*sqrt(576643*G^2 + 1805000*w^2)))/1729929)')
CalculatedLorWidth=(5000*(2673*w - sqrt(3)*sqrt(576643*G^2 + 1805000*w^2)))/1729929
disp('CalculatedGauWidth=sqrt(1729929*L^2 - 26730000*L*w + 25000000*w^2)/5000')
CalculatedGauWidth=sqrt(1729929*L^2 - 26730000*L*w + 25000000*w^2)/5000
disp(' ')
disp('If the alpha is unknown or varies from peak to peak, use shape 30'):
NumPeaks=1;
PeakShape=30; % Voigt with variable alpha
disp('          peak#       pos         height       width        area       alpha')
peakfit([x;v],0,0,NumPeaks,PeakShape,alpha)

% function v=voigt(x,pos,Gausswidth,voigtalpha)
% % Unit height Voigt profile function. x is the independent variable
% % (energy, wavelength, etc), Gausswidth is the Gaussian(Doppler) width,
% % and voigtalpha is ratio of the Gausswidth to the LorentzWidth (pressure
% % width). Version 3, August, 2019
% LorentzWidth=Gausswidth.*voigtalpha;
% if LorentzWidth<0, LorentzWidth=-LorentzWidth;end
% if Gausswidth<0, Gausswidth=-Gausswidth;end
% dx=x(2)-x(1); % x increment
% ex=[x-max(x)-dx x x+max(x)]; % Extended x
% gau=gaussian(ex,0,Gausswidth);
% lor=lorentzian(ex,pos,LorentzWidth);
% VoigtConv=ifft(fft(gau).*fft(lor))./sum(lor);
% g=VoigtConv./max(VoigtConv);
% oex=ex-max(x);
% outrange=val2ind(oex,0):val2ind(oex,max(x))-dx;
% v=g(outrange+1);