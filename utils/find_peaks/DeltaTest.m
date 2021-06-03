% Demonstration of the power spectrum of smoothing and differentation
% functions, obtained by applying them to a delta function.
smoothwidth=3;
smoothtype=1; % 1=rectangular; 2=triangular; 3=Gaussian
x=[1:100]';
delta=zeros(size(x));
delta(round(size(delta)./2))=1;delta(1)=0;
processeddelta=fastsmooth(delta,smoothwidth,smoothtype);
subplot(2,1,1)
plot(x,delta,'.')
xlabel('x');ylabel('y');
title('Delta Function')
subplot(2,1,2)
plot(x,processeddelta,'.')
title(['Processed Delta Function.   Smooth Width= ' num2str(smoothwidth) '   Smooth Type=' num2str(smoothtype) ])
xlabel('x');ylabel('y');