function g = BiGaussian(x,pos,wid,m)
% BiGaussian (different width on leading edge and trailing edge).
% pos=position; wid=width 
% m = ratio of widths of right-hand to left-hand halves.
% If m=1, it becomes identical to a Gaussian.
% Verison 2, T. C. O'Haver, 2012
% Example: Plots Gaussian and BiGaussian with m=3, both with halfwidth=20.
% x=[1:100];
% g=gaussian(x,50,20);
% bg=bigaussian(x,50,20,3);
% plot(x,g,x,bg)
%
lx=length(x);
hx=val2ind(x,pos);
hwid=2.*wid;
g(1:hx)=gaussian(x(1:hx),pos,hwid/(m+1));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,m.*hwid/(m+1));

function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);