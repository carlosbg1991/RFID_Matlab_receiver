% Script that demonstrates smoothing optimization for the measurement of
% the height of a noisy Gaussian peak.
% Requires fitgauss2.m, gaussfit.m, gaussian.m and fminsearch.m functions

x=0:.01:10;
SmoothType=1;
SmoothWidth=zeros(1,150);
SNR1=zeros(1,150);
SNR2=zeros(1,150);
SNR3=zeros(1,150);
    h1=zeros(1,150);h2=zeros(1,150);h3=zeros(1,150);
tic
for point=1:150,
    Width=1+2*point;
    y=exp(-(x-5).^2);
    NoiseArray=randn(1,1000);
    y=fastsmooth(y,Width,SmoothType,1);
    subplot(2,2,1)
    plot(x,y,'.')
    title('Gaussian, half width = 166 points')
    xlabel('SmoothType = 1')
    axis([0 10 0 1])
    SmoothedNoise=fastsmooth(NoiseArray,Width,3,SmoothType);
    PeakHeight=gaussfit(x(400:600),y(400:600));
    h1(point)=PeakHeight;
    h2(point)=std(SmoothedNoise);
    h3(point)=h1(point)./h2(point);
    SmoothWidth(point)=Width/166;
    
    subplot(2,2,2)
    plot(SmoothWidth,h1,'.')
    axis([0 2 0 1])
    xlabel('Smooth Ratio')
    ylabel('Peak Height')
    title('Peak Height')
    
    subplot(2,2,3)
    plot(SmoothWidth,h2,'.')
    axis([0 2 0 .5])
    ylabel('Standard deviation of the noise')
    xlabel('Smooth Ratio')
    title('Standard deviation of the noise')
    
    subplot(2,2,4)
    plot(SmoothWidth,h3,'.')
    axis([0 2 0 15])
    xlabel('Smooth Ratio')
    ylabel('Signal-to-Noise Ratio')
    title('Signal-to-Noise Ratio')
    drawnow
end
toc
