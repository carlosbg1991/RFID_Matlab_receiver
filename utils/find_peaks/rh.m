function h=rh(x,p,w)
% Reciprocal normalized Himmelblau function
% x is the independent variable, p is the position, and w is the width
% See https://en.wikipedia.org/wiki/Himmelblau%27s_function
% Example (for 10 different values of w):
% clf;hold on;x=-20:.1:20;for z=4:13;y=rh(x,0,z);plot(x,y);end;hold off
%
h=((x-p).^2+w-11).^2+((x-p)+w.^2-7).^2;
h=(min(h)./h);