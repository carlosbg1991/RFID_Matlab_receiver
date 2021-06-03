% Simple script using the peakfit.m function, showing how to
% extract the individual model components and the total model
% and plot them on a separate graph.
x=1:19;
y=[0 0 2 5 8 5 3 1 1 3 5 4 2 3 2 1 0 0 0];
% Fit a three Gaussian model to these data
[FitResults, GOF, baseline, coeff, residuals, xi, yi]=peakfit([x;y],0,0,3,1);
figure(2)
clf
% plot the individual model components, yi, and the total model, sum(yi). 
plot(xi,yi,xi,sum(yi))