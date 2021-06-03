% Script that shows that a Gaussian of unit height convoluted with a
% Gaussian of the same width is a Gaussian with a height of 1/sqrt(2) and
% a width of sqrt(2) and of equal area to the original Gaussian.
% Figure 2 shows an attempt to recover the original y from the convoluted
% yc by using the deconvgauss function. You can optionally add noise in
% line 9 to show how convolution smooths the noise and how deconvolution
% restores it. Requires gaussian.m, peakfit.m and deconvgauss.m in path.
% Version 2, March 2017
increment=.1;
Noise=.0;
width=100;
x=1:increment:700;
y=gaussian(x,351,width)+Noise.*randn(size(x)); % Original y
c=gaussian(x,351,width); % convolution function
yc=conv(y,c,'same')./sum(c); % yc= y convoluted with c
disp(' ')
disp('Demonstration of Gaussian convolution')

% Fit y and yc to a Gaussian function to show that it is a Gaussian and to
% determinie height, width and area
disp('Original y')
[FitResults,GOF]= peakfit([x;y]);
disp('           Peak      Position     Height       Width      Peak area');
disp(FitResults)
GoodnessOfFit=GOF

disp(' ')
disp('y convoluted with Gaussian')
[FitResults,GOF]= peakfit([x;yc]);
disp('          Peak      Position     Height       Width      Peak area');
disp(FitResults)
GoodnessOfFit=GOF

figure(1)
clf
plot(x,y,x,yc)
xlabel('x')
ylabel('y')
title('Blue: Original y     Green: yc = Convoluted y')

figure(2);
clf
% Attempt to recover the original y by deconvoluting a Gaussian of the
% same width
lx=length(x);
hlx=round(lx/2);
c=[gaussian(x(1:hlx),0,width) gaussian((x(hlx+1:lx)),max(x),width)];   % Gaussian convolution function, c
yr=ifft(fft(y)./fft(c)).*sum(c);
plot(x(400:6500),y(400:6500),x(400:6500),yr(400:6500))
xlabel('x')
ylabel('y')
title('Blue" Original y     Green: yr = recovered y')