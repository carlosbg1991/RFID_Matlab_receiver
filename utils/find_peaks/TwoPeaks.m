% Compare findpeaksG to peakfit.m with 2 very noisy peaks:
clf;
x=0:.01:10;
y=exp(-(x-3).^2)+exp(-(x-7).^2)+.2.*randn(size(x));
P=findpeaksG(x,y,1e-009,0.4,70,70,3);
sizeP=size(P);
NumPeaks=sizeP(1);
F=peakfit([x;y],0,0,NumPeaks,1);
disp('findpeaksG.m:');disp(P);
disp('peakfit.m:');disp(F);