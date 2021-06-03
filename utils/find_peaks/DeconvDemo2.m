% Deconvolution demo 2
xx=5:.1:200;
% Underlying signal with true (unknown) peak shape (Gaussian),heights, positions, and widths.
uyy=gaussian(xx,25,10);
% Observed signal yy, with noise added AFTER the broadening convolution (ExpG)
Noise=.00;
tc=300; % Time Constant
yy=expgaussian(xx,25,10,-tc)'+Noise.*randn(size(xx));
% Guess, or use prior knowledge, or curve fit one peak to ExpG to determine
% time constant (20), then compute transfer function cc
cc=exp(-(1:length(yy))./tc);   
% Attempt to recover original signal uyy by deconvoluting cc from yy
yydc=deconv([yy zeros(1,length(yy)-1)],cc).*sum(cc); % It's necessary to zero-pad the observed signal as shown here.
subplot(2,2,1);plot(xx,uyy);title('Underlying 4 Gaussian signal, uyy');
subplot(2,2,2);plot(xx,cc);title('Exponential transfer function, cc')
subplot(2,2,3);plot(xx,yy);title('observed broadened and noisy signal, yy');
subplot(2,2,4);plot(xx,yydc);title('After deconvoluting transfer function, yydc')
% Try more Noise (line 6) or a bad guess of time constant (line 11)
%  plot(xx,uyy,xx,yydc) % plot recovered signal overlaid with underlying signal
%  plot(xx,uyy,xx,yy) % plot observed signal overlaid with  with underlying signal
% Curve fit recovered signal to a Gaussian to determine peak parameters
% [FitResults,FitError]=peakfit([xx;yydc],26,42,1,1,0,10)
% Alternatively curve fit observed signal with an exponentially broadened
% Gaussian
% [FitResults,FitError]=peakfit([xx;yy],26,50,1,5,tc,10)