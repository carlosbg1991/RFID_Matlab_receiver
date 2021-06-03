function [Height, Position, Width]=lorentzfit(x,y)
% Converts y-axis to a reciprocal scale, fits a parabola
% (quadratic) to the (x,1/y) data, then calculates
% the position, width, and height of the
% Lorentzian from the three coefficients of the
% quadratic fit.  This is accurate only if the data have
% no baseline offset (that is, trends to zero far off the
% peak) and if there are no zeros or negative values in y.
% See gaussfit.m           (c) Tom O'Haver, 2013
%
% Example 1: Simplest Lorentzian data set
% [Height, Position, Width]=lorentzfit([1 2 3],[1 2 1]) 
%    returns Height = 2, Position = 2, Width = 2
%
% Example 2: best fit to synthetic noisy Lorentzian
% x=50:150;y=100.*lorentzian(x,100,100)+10.*randn(size(x));
% [Height,Position,Width]=lorentzfit(x,y) 
%   returns [Height,Position,Width] clustered around 100,100,100.
%
% Example 3: plots data set as points and best-fit Lorentzian as line
% x=[1 2 3 4 5];y=[1 2 3 2 1];
% [Height,Position,Width]=lorentzfit(x,y);
% plot(x,y,'o',linspace(0,8),Height.*lorentzian(linspace(0,8),Position,Width))
maxy=max(y);
for p=1:length(y)
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
z=ones(size(x))./y;
coef=polyfit(x,z,2);
Height=4*coef(1)./((4*coef(1)*coef(3))-coef(2)^2);
Position=-coef(2)/(2*coef(1));
Width=sqrt(((4*coef(1)*coef(3))-coef(2)^2)./coef(1))./sqrt(coef(1));

