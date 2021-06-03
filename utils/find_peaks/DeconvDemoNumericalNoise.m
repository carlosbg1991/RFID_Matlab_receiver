% Demonstration of the effect of numerical resolution on Gaussian
% convolution/deconvolution and the use of smoothing to reduce the
% numerical noise. 
% Create a rectangular function y, 400 points wide
increment=.01;
x=0:increment:20;
y=zeros(size(x));
y(800:1200)=1;           
y=rectanglepulse(x,10,5);

% Create a centered Gaussian convolution function (maximum at x=zero)
w=6; % w=width of the convolution function
c=gaussian(x,0,w)+gaussian(x,max(x),w);  % centered Gaussian convolution function

% Convolute rectangular function with Gaussian
yc=ifft(fft(y).*fft(c))./sum(c);  % Create Gaussian convoluted rectangular function

% Add a tiny bit of post-convolution random noise
yc=yc+.000000.*randn(size(yc)); % Noise added after the convolution

% Now attempt to recover the original signal by deconvolution (2 methods)
ydc=ifft(fft(yc)./fft(c)).*sum(c);  % Deconvolution by fft/ifft
% ydc=deconvgauss(x,yc,w);  % Deconvolution by "deconvgauss" function

% Plot all the steps
figure(1)
subplot(2,2,1); plot(x,y); title('original y');subplot(2,2,2); 
plot(x,c);title('Convolution function, c'); subplot(2,2,3); plot(x,yc(1:2001)); 
title('yc = y convoluted with c'); subplot(2,2,4); plot(x,ydc);title('ydc = recovered y, showing numerical noise');

% The numerical noise can be reduced by a little smoothing
SmoothWidth=8;
figure(2)
clf
plot(x,fastsmooth(ydc,SmoothWidth,1)) 
title(['Result of smoothing the recovered y with a Gaussian smooth of width ' num2str(SmoothWidth)] )