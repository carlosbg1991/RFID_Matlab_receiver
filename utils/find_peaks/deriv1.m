function d=deriv1(a)
% First derivative of vector using simple adjacent differences.
% Example: deriv([1 1 1 2 3 4]) returns [0 0 0 1 1 1]
%  T. C. O'Haver, 2013.
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
for j = 2:n;
  d(j)=a(j)-a(j-1);
end