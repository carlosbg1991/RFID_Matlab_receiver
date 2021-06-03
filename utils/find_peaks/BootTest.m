x=.001:.001:1;
Noise=.2.*randn(size(x));
SDofNoise=std(Noise)
y=1+x+Noise;
LeastSquaresMatlab
[Slope Intercept]
[SDSlope SDIntercept]
[coef, RSquared,BootResults]=plotit(x,y,1);
drawnow
pause
y=1+x+(y-1).*Noise;
LeastSquaresMatlab
[Slope Intercept]
[SDSlope SDIntercept]
[coef, RSquared,BootResults]=plotit(x,y,1);