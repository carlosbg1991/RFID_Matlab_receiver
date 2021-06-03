function d4=deriv4(a)
% 4th derivative of vector a
%  T. C. O'Haver, 2008.
% Example: Shows numerical round-off error in Matlab; a smaller increment
% than 0.0003 causes increasing noise. 
% x=-5:.0003:5; y=exp(-(x).^2);plot(deriv4(y))
n=length(a);
d4=zeros(size(a));
for j = 3:n-2;
  d4(j)=a(j+2) - 4.*a(j+1) + 6.*a(j) - 4.*a(j-1) + a(j-2);
end
d4(1:2)=d4(3);
d4(n-1:n)=d4(n-2);