function BroadenedPeak
% Matlab/Octave demo function comparing area of a Gaussian peak
% to the area of the same peak after exponential broadening.
Broadening=300; % time constant of exponential broadening function
x=-4:.005:11;
y=exp(-(x).^2); % Gaussian of unit height, peak at x=0; FWHM = 1.6651
y1=ExpBroaden(exp(-(x).^2)',-Broadening); % Apply broadening
OriginalArea=trapz(x,y)  % Trapezoidal numerical integration.
BroadenedArea=trapz(x,y1)  % Trapezoidal numerical integration.
figure(1)
plot(x,y,'.b',x,y1,'.r')
xlabel('Time')
ylabel('Detector signal')
title('Blue = Original peak      Red = broadened peak')
text(0,0.98.*max(y),['      <--- Area under curve = ' num2str(OriginalArea)])
text(0,0.99.*max(y1),['       <--- Area under curve = ' num2str(BroadenedArea)])
OriginalHeight=max(y)
BroadenedHeight=max(y1)
OriginalWidth=halfwidth(x,y)
BroadenedWidth=halfwidth(x,y1)
OriginalArea=OriginalArea
BroadenedArea=BroadenedArea
figure(2)
plot(x,cumsum(y),'.b',x,cumsum(y1),'.r')
xlabel('Time')
ylabel('Cumulative sum of detector signal')
title('Blue = Cumulative sum of original peak      Red = Cumulative sum of broadened peak')

function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant "t" by multiplying Fourier transforms and inverse
% transforming the result.
% Example:
% clg;
% x=1:100:
% y=zeros(size(x));
% y(40:60)=1;
% plot(x,y,x,ExpBroaden(y',-5))
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
