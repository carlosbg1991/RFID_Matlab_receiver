function d=deriv(a)
% Simple first derivative of vector using 2-point central difference.
%  T. C. O'Haver, 1988.
% Example:
%  x = [1 2 4 7 10 14 20 30];y = 2*x;plot(x,y,'o-');dy=deriv(x,y)
% ans =
%      2     3     5     6     7    10    16    20
% Note: to compute the derivative of y with respect to x, use derivxy.m
% x = [1 2 4 7 10 14 20 30];y = 2*x;derivxy(x,y)
% ans =
%      2     2     2     2     2     2     2     2
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end