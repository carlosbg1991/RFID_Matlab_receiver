function y=weibull(x,L,k)
% Weibull distributuion (normalized to height=1)
% See https://en.wikipedia.org/wiki/Weibull_distribution
% L is the peak x position. The FWHM of the peak is about 0.01737*k*L
% or k=FWHM/(.01737*L).
%
% Example 1:
% x=0:.01:4.5; 
% k=15;
% y=weibull(x,.8,12)+weibull(x,2,30)+weibull(x,3.5,50);
% plot(x,y);
%
% Example 2:
% x=0:.01:2.5;
% L=1;
% clf;
% hold on;
% for k=1:.2:5;
%     plot(x,weibull(x,L,k));
% end;
% hold off
%
if x<0,
    y=0;
else
y=(k/L).*(x/L).^(k-1).*exp(-(x/L).^k);
y=y./max(y);
end
