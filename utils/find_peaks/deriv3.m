function d3=deriv3(a)
% Third derivative of vector a
%  T. C. O'Haver, 2008.
n=length(a);
d3=zeros(size(a));
for j = 3:n-2;
  d3(j)=a(j+2) - 2.*a(j+1) + 2.*a(j-1) - a(j-2);
end
d3(1:2)=d3(3);
d3(n-1:n)=d3(n-2);
