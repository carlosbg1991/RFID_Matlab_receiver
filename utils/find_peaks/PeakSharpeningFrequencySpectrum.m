x=[0:.01:18];
y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2);
Enhancedsignal=enhance(y,1000,1000000,3);
subplot(2,1,1);
plot(x,y,x,Enhancedsignal,'r')
title('Four Gaussian peaks of equal width')
xlabel('Blue: Original signal     Red: After peak sharpening, showing improved resolution.')
ylabel('y')
fy=fft(y);
sy=fy .* conj(fy);
plotrange=1:length(fy)/2;
realsy=real(sy(plotrange));
f=((plotrange-1)./range(x));
fs=([f;realsy]);
subplot(2,1,2);
semilogx(fs(1,:),fs(2,:))

hold on

fy=fft(Enhancedsignal);
sy=fy .* conj(fy);
plotrange=1:length(fy)/2;
realsy=real(sy(plotrange));
f=((plotrange-1)./range(x));
fs=([f;realsy]);

semilogx(fs(1,:),fs(2,:),'r-')
title('Frequency Spectra of the above signals')
xlabel('Blue: Original signal     Red: After peak sharpening, showing emphasis of higher frequencies.')
ylabel('semilogx plot')
hold off
  