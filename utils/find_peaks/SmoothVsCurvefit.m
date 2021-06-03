% Comparison of measuring peak height by smoothing and by curve fitting.
% The peak height, position and width of the peak is specified by line 7-8.
% The smooth type is specified in  line 9: 1=rectangle; 2=triangle; 3=pseudo
% Gaussian. Requires peakfit.m and gaussian.m in path.
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
    ny=y+.1*randn(size(x));
    fitresults=peakfit([x;ny]);
    peakfitheight(trial)=fitresults(3);
     plot(x,fastsmooth(ny,SmoothWidth(trial),SmoothType,1))
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
clg
plot(SmoothRatio,Maxy,SmoothRatio,SmoothedNoise,SmoothRatio,SNR,SmoothRatio,peakfitheight)
switch SmoothType
    case 1
        title(' Comparison of measuring peak height by curve fitting and by smoothing with Rectangular (boxcar) smooth')
    case 2
        title('Comparison of measuring peak height by curve fitting and by smoothing with Triangular smooth')
    case 3
        title('Comparison of measuring peak height by curve fitting and by smoothing with Pseudo-Gaussian smooth')
end
xlabel('Smooth Ratio = SmoothWidth/PeakWidth')
ylabel('Blue/Cyan = peak heights.    Green = Noise.    Red = Signal-to-noise ratio.')
text(max(SmoothRatio)-1,SNR(trial),'Signal-to-noise ratio (SNR)')
text(max(SmoothRatio)-1,Maxy(trial),'Peak height after smoothing')
text(max(SmoothRatio)-1,SmoothedNoise(trial),'Noise remaining in smooth signal')
text(max(SmoothRatio)-1,peakfitheight(trial),'peak height by curve fitting unsmoothed signals')