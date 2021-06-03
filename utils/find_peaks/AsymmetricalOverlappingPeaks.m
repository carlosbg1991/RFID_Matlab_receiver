% Mystery Peak script
% A combination of first-derivative symmetrization and curve fitting is
% used to analyze a complex mystery peak.
% Requires peakfit.m, symmetrize.m, halfwidth.m, val2ind.m, and
% fastsmooth.m in the path. Download from
% https://terpconnect.umd.edu/~toh/spectrum/
%
% To use a dataset identical to the one on the plots:
load MysteryPeak 
% or download and run MysteryPeak.m to generate another dataset with the
% same peak but independent noise samples:
% [x,y]=MysteryPeak;

disp('--------------------------------------------------------------------')

disp(' ')
figure(1)
disp('Figure 1: Two exponentially broadened Gaussians (variable shape #39)')
[FitResults,GOF]=peakfit([x y],0,0,2,39,1);
disp('          peak#         pos         height     width        area       lambda')
disp(FitResults)
disp(['Percent fitting error: ' num2str(GOF(1)) '     R2: '  num2str(GOF(2)) ])
subplot(2,1,1);title('Two exponentially broadened Gaussians (variable shape #39)')
disp(' ')

figure(2)
disp('Figure 2: Three Gaussians model')
[FitResults,GOF]=peakfit([x y],0,0,3,1);
disp('          peak#         pos         height     width        area')
disp(FitResults)
disp(['Percent fitting error: ' num2str(GOF(1)) '     R2: '  num2str(GOF(2)) ])
subplot(2,1,1);title('Three Gaussians model')
disp(' ')

figure(3)
disp('Figure 3: Four Gaussians model ')
[FitResults,GOF]=peakfit([x,y],0,0,4,1,0,10);
disp('          peak#         pos         height     width        area   ')
disp(FitResults)
disp(['Percent fitting error: ' num2str(GOF(1)) '     R2: '  num2str(GOF(2)) ])
subplot(2,1,1);title('Four Gaussians model');
disp(' ')

figure(4)
disp('Figure 4: Determining the required smooth width from the number of')
disp('points in the halfwidth of the original peak and using one-tenth of that.')
clf;
WidthPoints=halfwidth(1:length(x),y); % Number of data points in halfwidth of peak.
% Smooth ratio of 0.1 produces no distortion
SmoothWidth=0.1.*WidthPoints
smy=fastsmooth(y,SmoothWidth,3);
plot(x,derivxy(x,y),'g',x,y,'k',x,smy,'c',x,derivxy(x,smy),'r')
title('Black: original.  Cyan: smoothed.  Green:  1st derivative.  Red: smoothed 1st derivative')
xlabel('Based on the amplitude of the derivative, the optimum value of symmfactor is near 1')

disp(' ')
figure(5)
disp('Figure 5: Determining the optimum symmetrization factor (symmfactor)')
clf;
hold on;
inx1=val2ind(x,58); % start index of zoom-in range for trailing edge
inx2=val2ind(x,65); % end index of zoom-in range for trailing edge
for symmfactor=1:.2:2
    sy=symmetrize(x,y,symmfactor,SmoothWidth,3,1)';
    plot(x(inx1:inx2),sy(inx1:inx2))
end
hold off
grid
title('Trailing edge. SymmFactors from 1 (blue) to 2 (cyan), in steps of 0.1.')
xlabel('This graph shows that the optimum value of symmfactor is 1.3 (yellow line)')
ylabel('Smoothed symmetrized signal')

disp(' ')
figure(6)
disp('Figure 6: Symmetrizing the peak and fitting to a 3-Gaussian model')
symmfactor=1.25; % From analysis of Figure 5
disp('Three Gaussians model after symmetrization')
tic
symy=symmetrize(x,y,symmfactor,SmoothWidth,3,1)';
[FitResults,GOF]=peakfit([x symy],57,15.5,3,1,0,10);
toc
disp('          peak#         pos         height     width        area')
disp(FitResults)
disp(['Percent fitting error: ' num2str(GOF(1)) '     R2: '  num2str(GOF(2)) ])
subplot(2,1,1);
title(['Three Gaussians model after symmetrization.  SmoothWidth: ' num2str(SmoothWidth) '   SymmFactor: ' num2str(num2str(symmfactor)) ])
disp(' ')

figure(7)
disp('Figure 7: Direct fitting of a three equal-alpha exponential ')
disp('Gaussian model to the original data gives less accurate results.')
alpha=symmfactor./(x(2)-x(1)) % alpha is the time constant expressed in number of data points
tic
[FitResults,GOF]=peakfit([x y],58,30,3,5,alpha);
toc
disp('          peak#         pos         height     width        area')
disp(FitResults)
disp(['Percent fitting error: ' num2str(GOF(1)) '     R2: '  num2str(GOF(2)) ])
subplot(2,1,1);title('Three Exponential Gaussians model applied to original data')
disp(' ')