% Demonstration of fitting a set of simulated spectra with fixed-position
% Gaussian peaks, with variable height.
FitResults=zeros(100,2,5);
y=zeros(100,101);
x=0:.1:10;
for time=1:100,
    y(time,:)=time.*exp(-(x-5).^2)+50*exp(-(x-3).^2)+randn(size(x));
end;
 for time=1:100;
     % Shape 16 is fixed-position Gaussian; positions are specified in [3 5]
     FitResults(time,:,:)=peakfit([x;y(time,:)],0,0,2,16,0,0,0,0,[3 5]);
     drawnow;
 end
 clf
plot(1:time,FitResults(1:time,2,3),1:time,FitResults(1:time,1,3))
xlabel('Time')
ylabel('Position')
text(0,50,'   Peak 1')
text(0,5,'   Peak 2')
axis([1 100 0 100])
title ('Peak Heights of two fixed-position Gaussians (shape 16)')