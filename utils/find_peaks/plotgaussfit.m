function [Height, Position, Width]=plotgaussfit(x,y)
% Converts y-axis to a log scale, fits a parabola (quadratic) to the
% (x,ln(y)) data, calculates the position, width, and height of the
% Gaussian from the three coefficients of the quadratic fit, then plots the
% raw data as red dots and the best-fit Gaussian as a line.  This works
% only if the data have no baseline offset (that is, trends to zero far off
% the peak) and if there are no zeros or negative values in y. Example:
% [Height, Position, Width]=plotgaussfit([1 2 3],[1 2 1]) returns Height = 2,
% Position = 2, Width = 2
maxy=max(y);
for p=1:length(y),
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
z=log(y);
coef=polyfit(x,z,2);
a=coef(3);
b=coef(2);
c=coef(1);
Height=exp(a-c*(b/(2*c))^2);
Position=-b/(2*c);
Width=2.35482/(sqrt(2)*sqrt(-c));
xx=linspace(min(x), max(x));
plot(x,y,'r.',xx,Height.*gaussian(xx,Position,Width))

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);