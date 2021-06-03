% Demonstration of the effect of smoothing on the peak height, random 
% white noise, and signal-to-noise ratio of a peak specified by line 8.
% The smooth type is specified by line 9: 1=rectangle; 2=triangle; 
% 3=pseudo Gaussian
interval=.5;
x=1:interval:200;
PeakWidth=20;
y=gaussian(x,100,PeakWidth); % You may substitute Lorentzian or other peak function
SmoothType=3; %  1=rectangle; 2=triangle; 3=pseudo Gaussian
NoiseArray=randn(1,10000); % Long array here gives more consistent results.
for trial=1:round(50/interval);
    SmoothWidth(trial)=round(trial);
    smoothedy=fastsmooth(y,SmoothWidth(trial),SmoothType);
    Maxy(trial)=max(smoothedy);
    SmoothedNoise(trial)=std(fastsmooth(NoiseArray,SmoothWidth(trial),SmoothType));
    SNR(trial)=Maxy(trial)./SmoothedNoise(trial);
    SmoothRatio(trial)=SmoothWidth(trial)./(PeakWidth./interval);
     plot(x,fastsmooth(y+.1*randn(size(x)),SmoothWidth(trial),SmoothType,1))
     axis([0 200 -.1 1]);
     drawnow
     xlabel('x')
     ylabel('Smoothed y')
     title([ ' Smooth Ratio = ' num2str(SmoothRatio(trial)) ] )
     if trial==1;
         title('Original signal before smoothing')
         pause(2)
     else
         pause(.05)
     end
end
plot(SmoothRatio,Maxy,SmoothRatio,SmoothedNoise,SmoothRatio,SNR)
switch SmoothType
    case 1
        title('Rectangular (boxcar) smooth')
    case 2
        title('Triangular smooth')
    case 3
        title('Pseudo-Gaussian smooth')
end
xlabel('Smooth Ratio = SmoothWidth/PeakWidth')
ylabel('Blue = peak height.    Green = Noise.    Red = Signal-to-noise ratio.')
text(max(SmoothRatio)-.6,SNR(trial),'Signal-to-noise ratio')
text(max(SmoothRatio)-.5,Maxy(trial),'Peak height')
text(max(SmoothRatio)-.5,SmoothedNoise(trial),'Noise')