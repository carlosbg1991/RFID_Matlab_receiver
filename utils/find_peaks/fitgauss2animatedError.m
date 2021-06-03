function err = fitgauss2animatedError(lambda,t,y)
% Fitting function for multiple overlapping Gaussians, with statements
% added (lines 18 and 19) to slow the progress and plot each step along the
% way, for educational purposes.
%  T. C. O'Haver, 2006
global c NumTrials TrialError
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end
c = A\y';
z = A*c;
err = norm(z-y'); 
peak1=c(1).*gaussian(t,lambda(1),lambda(2));
peak2=c(2).*gaussian(t,lambda(3),lambda(4));
model=peak1+peak2;
subplot(1,2,1)
plot(t,y,'ro',t,peak1,'c',t,peak2,'m',t,model,'b')
title(num2str(err))
NumTrials=NumTrials+1;
TrialError(NumTrials)=err;
subplot(1,2,2)
semilogy(1:NumTrials,TrialError)
drawnow
% pause(.01)