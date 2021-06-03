function err = fitgauss2b(lambda,t,y)
% Fitting functions for a Gaussian band spectrum with flat 
% baseline correction. Returns baseline as c(1)
%  T. C. O'Haver, 2013   
global c
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end
B=[ones(size(y))' A];
c = B\y';
z = B*c;
err = norm(z-y');