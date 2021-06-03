% Demonstration of Gaussian convolution and deconvolution of Gaussian peak.
% Requires Gaussian and fastsmooth functions .

increment=.01;
position=10;
width=.5; % Width of the original Gaussian peak.
cw=2; % cw = Convolution width of the physical convolution process (unknown)
dw=2.02; % dw = Deconvolution width (estimated) must equal cw for perfect results
SmoothWidth=10; % Adjust for best results to reduce noise.
%
% Create simulated signal
x=0:increment:20;
y=gaussian(x,position,width);  % Create a Gaussian function y of width w.        

% Add pre-convolution random noise
y=y+.01.*randn(size(y));    % Noise added before the convolution
c=gaussian(x,0,cw)+gaussian(x,max(x),cw);  % zero centered Gaussian convolution function
widths=[ 1.97 1.98 1.99 2.00 2.01 2.02 2.03];
for trial=1:7,
% Create a zero centered Gaussian convolution function (maximum at x=zero)
dw=widths(trial); % dw = Deconvolution width (estimated) must equal cw for perfect results

c2=gaussian(x,0,dw)+gaussian(x,max(x),dw);  % zero centered Gaussian deconvolution function

% Convolute rectangular function with Gaussian
yc=ifft(fft(y).*fft(c))./sum(c);  % Create Gaussian convoluted rectangular function

% Add a little bit of post-convolution random white noise
yc=yc+.000001.*randn(size(yc)); % Noise added after the convolution

% Now attempt to recover the original signal by deconvolution 
ydc=ifft(fft(yc)./fft(c2)).*sum(c2);  % Deconvolution by fft/ifft

sy=fastsmooth(ydc,SmoothWidth,3); % Smooth the result with smooth width=SmoothWidth
sy=fastsmooth(sy,SmoothWidth,3); % Smooth again 

% Plot all the steps
subplot(2,2,1); plot(x,y); title('Original signal,y');
xlabel(['FWHM = ' num2str(halfwidth(x,y,position)) ]);
subplot(2,2,2); plot(x,yc(1:2001)); title(['yc=y convoluted with Gaussian of width = ' num2str(cw) ]); 
xlabel(['FWHM = ' num2str(halfwidth(x,yc,position)) ]);
subplot(2,2,3); plot(x,ydc);title(['ydc=Recovered y. Deconvolution width = ' num2str(dw) ]);
subplot(2,2,4); plot(x,sy);title('Smoothed recovered y');
xlabel(['FWHM = ' num2str(halfwidth(x,sy,position)) ]);
pause(1)
end