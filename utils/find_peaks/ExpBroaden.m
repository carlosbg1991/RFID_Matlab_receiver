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