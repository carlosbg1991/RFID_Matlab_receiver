% Demonstration of Gaussian convolution and deconvolution of two
% closely-spaced underlying lorentzian peaks. Requires lorentzian,
% and fastsmooth functions .

increment=.01;
position=10;
width=.5; % Width of the original lorentzian peak.
cw=1.1; % cw = Convolution width of the physical convolution process (unknown)
dw=1.1; % dw = Deconvolution width (estimated) must equal cw for perfect results
SmoothWidth=5; % Adjust for best results to reduce noise.
%
% Create simulated signal with two closely-spaced underlying peaks, height
% ratio 2:1
x=0:increment:20;
y=1.3.*lorentzian(x,position,width)+1.3./2.*lorentzian(x,position+1,width);  % Create a lorentzian function y of width w.        

% Add pre-convolution random noise
y=y+.01.*randn(size(y));    % Noise added BEFORE the convolution

% Create a zero centered Gaussian convolution function (maximum at x=zero, width=cw)
c=gaussian(x,0,cw)+gaussian(x,max(x),cw);  % zero centered Gaussian convolution function
% Create a zero centered Gaussian deconvolution function (maximum at x=zero, width=c2)
c2=gaussian(x,0,dw)+gaussian(x,max(x),dw);  % zero centered Gaussian deconvolution function

% Convolute underlying peaks with Gaussian convolution function
yc=ifft(fft(y).*fft(c))./sum(c);  % Create Gaussian convoluted rectangular function

% Add a little bit of post-convolution random white noise to create
% simulated observed signal
yc=yc+.000001.*randn(size(yc)); % Noise added AFTER the convolution

% Now attempt to recover the original signal by deconvolution 
ydc=ifft(fft(yc)./fft(c2)).*sum(c2);  % Deconvolution by fft/ifft

sy=fastsmooth(ydc,SmoothWidth,3); % Smooth the result with smooth width=SmoothWidth
sy=fastsmooth(sy,SmoothWidth,3); % Smooth again 

% Plot all the steps
figure(1)
subplot(2,2,1); plot(x,y); title('Underlying peaks, y, not observable');
xlabel(['FWHM = ' num2str(halfwidth(x,y,position)) ]);
subplot(2,2,2); plot(x,yc(1:2001)); title(['Observed signal yc after convolution with Gaussian of width = ' num2str(cw) ]); 
xlabel(['FWHM = ' num2str(halfwidth(x,yc,position)) ]);
subplot(2,2,3); plot(x,ydc);title(['Raw recovered signal, ydc, after deconvolution width = ' num2str(dw) ]);
subplot(2,2,4); plot(x,sy);title('Smoothed ydc');
xlabel(['FWHM = ' num2str(halfwidth(x,sy,position)) ]);
