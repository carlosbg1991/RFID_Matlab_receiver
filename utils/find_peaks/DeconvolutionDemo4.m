% Deconvolution demo with 4 Gaussian peaks.
% Similar to the first figure on
% https://terpconnect.umd.edu/~toh/spectrum/Deconvolution.html
xx=0:.1:255;
% Underlying signal with true peak shape (Gaussian),heights, positions, and widths.
uyy=1.2.*gaussian(xx,25,10)+gaussian(xx,55,10)+0.8.*gaussian(xx,85,10)+0.7.*gaussian(xx,115,10);   
% Observed signal yy,
Noise=.0; % amount of noise added BEFORE the broadening convolution
Noise2=.0001; % amount of noise added AFTER the broadening convolution
uyy=uyy+Noise.*randn(size(xx));
tc=120; % Exponential facotr (Time constant)
yy=expbroaden(uyy',-tc)'; % Convolution with exponential function
yy=yy+Noise2.*randn(size(xx));
% Guess, or use prior knowledge, or curve fit one peak to ExpG to determine
% time constant (20), then compute transfer function cc
cc=exp(-(1:length(yy))./tc);   
% Attempt to recover original signal uyy by deconvoluting cc from yy
yydc=deconv([yy zeros(1,length(yy)-1)],cc).*sum(cc); % It's necessary to zero-pad the observed signal as shown here.
subplot(2,2,1);plot(xx,uyy);title('Underlying Gaussian signal, uyy');
subplot(2,2,2);plot(xx,cc);title('Exponential transfer function, cc')
subplot(2,2,3);plot(xx,yy);title('observed broadened and noisy signal, yy');
subplot(2,2,4);plot(xx,yydc);title('After deconvoluting transfer function, yydc')

% Use the findpeaksG function to measure peak position, heights, widths
disp('Peaks of the Underlying Gaussian signal, uyy')
disp('          Peak#     Position     Height      Width          Area')
Puyy=findpeaksG(xx,uyy,6.1467e-005,0.6,10,10,3);
disp(Puyy)
disp(' ')
disp('Peaks of observed broadened and noisy signal, yy')
disp('          Peak#     Position     Height      Width          Area')
Pyy=findpeaksG(xx,yy,6e-005,0.2,10,10,3);
disp(Pyy)
disp(' ')
disp('Peaks of the recovered signal, yydc')
disp('          Peak#     Position     Height      Width          Area')
Pyydc=findpeaksG(xx,yydc,6.1467e-005,0.6,50,80,3);
disp(Pyydc)
 