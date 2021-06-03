% Script to compare shapes of derivatives for a Gaussian peak. Computes and
% plots the normalized amplitudes of the zeroth through fourth derivatives.
% Requires deriv.m and gaussian.m funcitons in path.
% See https://terpconnect.umd.edu/~toh/spectrum/functions.html for other
% shape functions other than Gaussian to use in line 8. 
clf
x=-5:100;
y=gaussian(x,50,25);
d1=deriv(y);
d1n=(d1-min(d1))./max(d1-min(d1));
d2=deriv(d1);
d2n=(d2-min(d2))./max(d2-min(d2));
d3=deriv(d2);
d3n=(d3-min(d3))./max(d3-min(d3));
d4=deriv(d3);
d4n=(d4-min(d4))./max(d4-min(d4));
plot(x,y,x,d1n,x,d2n,x,d3n,x,d4n)
axis([0 100 0 1])
title('Derivative orders: 0=blue     1=green     2=red      3=cyan     4-magenta')
xlabel('x')
ylabel('Normalized numerical amplitude')
