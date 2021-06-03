function err = fitgauss(lambda,t,y)
%   Fitting functions for multiple overlapping Gaussians
%  T. C. O'Haver, 2006  
global c
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end
c = A\y';
z = A*c;
err = norm(z-y');