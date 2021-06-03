function d=derivxy(x,y)
% First derivative of y with respect to x using 2-point central difference.
%  T. C. O'Haver, 2011. Corrected 2/4/14
% Example:
%  x = [1 2 4 7 10 14 20 30];y = 2*x;plot(x,y,'o-');dy=derivxy(x,y)
% ans =
%      2     2     2     2     2     2     2     2
% Compare to:
% x = [1 2 4 7 10 14 20 30];y = 2*x;deriv(y)
% ans =
%      2     3     5     6     7    10    16    20
%
n=length(y);
d=zeros(size(y));
d(1)=(y(2)-y(1))./(x(2)-x(1));
d(n)=(y(n)-y(n-1))./(x(n)-x(n-1));
for j = 2:n-1
  d(j)=(y(j+1)-y(j-1)) ./ ((x(j+1)-x(j-1)));
end
