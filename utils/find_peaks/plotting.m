% Examples of creating and plotting data sets in matrix form
%
% FIRST, you must download the "gaussian" function from 
% http://terpconnect.umd.edu/~toh/spectrum/gaussian.m
% and install it on your Matlab/Octave path.

clear A
clf

w=200:800;
a=gaussian(w,500,100);
subplot(2,2,1);
plot(w,a)
title('Plotting a single Gaussian curve')

for n=1:10;
    A(:,n)=n.*gaussian(w,500,100);
end;
subplot(2,2,2);
plot(w,A)
title('Gaussians with 10 different heights')

w=200:800;
for n=1:10;
    A(:,n)=gaussian(w,200+(50*n),100);
end;
subplot(2,2,3);
plot(w,A)
title('Gaussians with 10 different positions')

w=200:800;
for n=1:10;
    A(:,n)=gaussian(w,500,20+(20*n));
end;
subplot(2,2,4);
plot(w,A)
title('Gaussians with 10 different widths')