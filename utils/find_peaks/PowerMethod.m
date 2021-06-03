% Demonstration of power method suggested by M. Farooq Wahab
xx=0:.1:128;
power=5;
SmoothWidth=5; % Smooth width for convoluted signal (e.g. 3 - 12)
% Underlying signal with true peak shape (Gaussian),heights, positions, and widths.
uyy=1.2.*gaussian(xx,25,10)+gaussian(xx,55,10);   
% Observed signal yy, with  (ExpG)
Noise=.001; % amount of noise added BEFORE the broadening convolution
Noise2=.001; % amount of noise added AFTER the broadening convolution
tc=120; % Exponential facotr (Time constant)
%
uyy=uyy+Noise.*randn(size(xx));
yy=expbroaden(uyy',-tc)'; % Convolution with exponential function
yy=yy+Noise2.*randn(size(xx));
% Guess, or use prior knowledge, or curve fit one peak to ExpG to determine
% time constant (20), then compute transfer function cc
cc=exp(-(1:length(yy))./tc);  
yyp=yy.^power;
% Attempt to recover original signal uyy by deconvoluting cc from yy
figure(1)
yydc=deconv([yy zeros(1,length(yy)-1)],cc).*sum(cc); % It's necessary to zero-pad the observed signal as shown here.
% Smooth to reduce blue noise amplified by deconvolution
yydc=fastsmooth(yydc,SmoothWidth,3);
subplot(2,2,1);plot(xx,uyy);title('Underlying Gaussian signal, uyy');
subplot(2,2,2);plot(xx,cc);title('Exponential transfer function, cc')
subplot(2,2,3);plot(xx,yy);title('observed broadened and noisy signal, yy');
subplot(2,2,4);plot(xx,yyp);title(['Raising the signal to a power, yyp. Power= ' num2str(power)] )

% Use the findpeaksG function to measure peak position, heights, widths
disp('Peaks of the Underlying Gaussian signal, uyy')
disp('          Peak#     Position     Height      Width          Area')
Puyy=findpeaksG(xx,uyy,6e-005,0.6,10,10,3);
disp(Puyy)
SeparationToWidthRatio=(Puyy(2,2)-Puyy(1,2))./Puyy(1,4)
disp(' ')
disp('Peaks of observed broadened and noisy signal, yy')
disp('          Peak#     Position     Height      Width          Area')
Pyy=findpeaksG(xx,yy,6e-005,0.2,10,10,3);
disp(Pyy)
SeparationToWidthRatio=(Pyy(2,2)-Pyy(1,2))./Pyy(1,4)
disp(' ')
disp('Peaks of the smoothed deconvoluted recovered signal, yydc')
disp('          Peak#     Position     Height      Width          Area')
Pyydc=findpeaksG(xx,yydc,6e-005,0.6,50,80,3);
disp(Pyydc)
SeparationToWidthRatio=(Pyydc(2,2)-Pyydc(1,2))./Pyydc(1,4)
 disp(' ')
disp('Peaks of the observed broadened and noisy signal raised to power, yyp')
disp('          Peak#     Position     Height      Width          Area')
Pyy2=findpeaksG(xx,yyp,1e-005,0.01,10,10,3);
disp(Pyy2)
SeparationToWidthRatio=(Pyy2(2,2)-Pyy2(1,2))./Pyy2(1,4)
figure(2)
clf
plot(xx,uyy./max(uyy),xx,yydc./max(yydc),xx,yyp./max(yyp))
title('Blue: Original signal; Green: Deconvolution method; Red: Power method, both normalized to 1.00')
xlabel('x')
ylabel('Normalized signals')
text(1,1.15,['Power method: n=' num2str(power)])
text(1,1.1,['Deconvolution: SmoothWidth=' num2str(SmoothWidth)])