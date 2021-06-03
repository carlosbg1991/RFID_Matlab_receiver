format compact
format short g
deltax=.01;
x=-5:deltax:5;
mu=0;
sigma=1;
disp(' ')
disp('ANALYTICAL DERIVATIVES')
% Analytical expressions for Gaussian and its derivatives (from Wolfram
% Alpha, e.g. http://www.wolframalpha.com/input/?i=second+derivative+of+gaussian)
% g=Gaussian
g= exp(-(x - mu).^2./(2.*sigma.^2))./(sqrt(2.*pi).*sigma);
% gd2=2nd derivative of Gaussian
gd2=(((x-mu).^2.*exp(-(x-mu).^2./(2*sigma.^2)))./sigma.^4-exp(-(x - mu).^2./(2*sigma.^2))./sigma.^2)./(sqrt(2*pi)*sigma);
% gd4=4th derivative of Gaussian
gd4=((((x - mu).^4.*exp(-(x - mu).^2./(2*sigma.^2)))./sigma.^8 - (6*(x - mu).^2.*exp(-(x - mu).^2./(2.*sigma.^2)))./sigma.^6 + (3.*exp(-(x - mu).^2./(2.*sigma.^2)))./sigma.^4)./(sqrt(2.*pi)*sigma));
figure(1);clf
plot(x,g,x,gd2,x,gd4)
title('Theoretical analytical Gaussian peak (blue) and its 2nd (red) and 4th (yellow) derivatives')
xlabel('X');ylabel('Amplitude')

% Numerical derivatives
disp(' ')
disp('NUMERICAL DERIVATIVES by repeated application of first differentation')
% Second derivative is two successive first differentations of the original
% Gaussian g
nd2=derivxy(x,derivxy(x,g));
figure(2)
plot(x,gd2,x,nd2,'.')
title('Numerical 2 x 1st derivatives (red dots) vs analytical 2nd derivative (Blue line) ')
xlabel('X');ylabel('Amplitude')
PercentDifference2nd=100*(max(gd2)-max(nd2))/max(nd2)
figure(3)
% Fourth derivative is two successive first differentations of the second derivative ndg.
nd4=derivxy(x,derivxy(x,nd2));
plot(x,gd4,x,nd4,'.');
title('Numerical 2 x 2nd derivatives (red dots) vs analytical 4th derivative (Blue line)')
xlabel('X');ylabel('Amplitude')
PercentDifference4th=100*(max(gd4)-max(nd4))/max(nd4)
figure(4)
 plot(x,gd4,x,deriv4(g)/deltax^4,'.')
title('Numerical 4nd derivative in one step using deriv4.m (red dots) vs analytical 4th derivative (Blue line)')
xlabel('X');ylabel('Amplitude')