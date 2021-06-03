% Simple script for a Monte Carlo simulation of the measurement of the
% slope of a noisy x,y straight-line data set with slope=1, intercept=0, and
% added random n oise with a standard deviaiton of 5& oif the maximum y
% value
x=1:50;
for trial=1:50;
    y=1+x+5.*randn(size(x));
    [coef, RSquared]=plotit(x,y,1);
    slopes(trial)=coef(1);
    axis([-5 55 -5 55]);
    drawnow;
    pause(.1)
end
AverageSlope=mean(slopes)
StandardDeviation=stdev(slopes)
