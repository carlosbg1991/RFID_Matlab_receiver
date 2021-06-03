% Simple Matlab script for calculating the first-order least-square fit of y vs x,
% including the Slope and Intercept and the predicted standard deviation of
% the slope (SDSlope) and intercept (SDIntercept).

NumPoints=length(x);
Sxx = sum((x-mean(x)).^2);
Syy = sum((y-mean(y)).^2);
Sxy = sum((x-mean(x)).*(y-mean(y)));
Slope = Sxy./Sxx;
Intercept = mean(y)-Slope*mean(x);
Sy = sqrt((Syy-Slope^2*Sxx)/(NumPoints-2));
SDSlope = Sy/sqrt(Sxx);
SDIntercept = Sy*sqrt(1./(NumPoints-(sum(x).^2)./sum(x.^2)));

% Optional plot of data and fitted line
plotx = linspace(min(x),max(x));
plot(x,y,'o',plotx,Slope*plotx +Intercept)