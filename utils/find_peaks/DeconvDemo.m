% Deconvolution demo
xx=5:.1:65;
% Underlying signal with true (unknown) peak shape (Gaussian),heights, positions, and widths.
uyy=modelpeaks2(xx,[1 1 1 1],[1.2 1.1 1 .9],[10 20 30 40],[3 4 5 6],[0 0 0 0]);
% Observed signal yy, with noise added AFTER the broadening convolution (ExpG)
Noise=.001;
yy=modelpeaks2(xx,[5 5 5 5],[1.2 1.1 1 .9],[10 20 30 40],[3 4 5 6],[-20 -20 -20 -20])+Noise.*randn(size(xx));
% Guess, or use prior knowledge, or curve fit one peak to ExpG to determine
% time constant (20), then compute transfer function cc
cc=exp(-(1:length(yy))./20);   
% Attempt to recover original signal uyy by deconvoluting cc from yy
yydc=deconv([yy zeros(1,length(yy)-1)],cc).*sum(cc); % It's necessary to zero-pad the observed signal as shown here.
subplot(2,2,1);plot(xx,uyy);title('Underlying 4 Gaussian signal, uyy');
subplot(2,2,2);plot(xx,cc);title('Exponential transfer function, cc')
subplot(2,2,3);plot(xx,yy);title('observed broadened and noisy signal, yy');
subplot(2,2,4);plot(xx,yydc);title('After deconvoluting transfer function, yydc')
% Try more Noise (line 6) or a bad guess of time constant (line 10)
%  plot(xx,uyy,xx,yydc) % plot recovered signal overlaid with underlying signal
%  plot(xx,uyy,xx,yy) % plot observed signal overlaid with  with underlying signal
% Curve fit recovered signal to 4 Gaussians to determine peak parameters
% [FitResults,FitError]=peakfit([xx;yydc],26,42,4,1,0,10)
% Alternatively curve fit observed signal with 4 exponentially broadened Gaussians
% [FitResults,FitError]=peakfit([xx;yy],26,50,4,5,20,10)