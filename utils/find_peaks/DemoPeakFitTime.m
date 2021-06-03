% Demonstration of fitting a set of simulated spectra with 
% Gaussian peaks with variable positions.
FitResults=zeros(100,2,5);
y=zeros(100,101);
x=0:.1:10;
for time=1:100,
    y(time,:)=exp(-(x-(5+time/30)).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
end;
 for time=1:100;
     Results=peakfit([x;y(time,:)],0,0,2,1,0,2,[3 1 6 1]);
     FitResults(time,:,:)=Results;
     drawnow;
 end
 clf
plot(1:time,FitResults(1:time,2,2),1:time,FitResults(1:time,1,2),1:time,FitResults(1:time,2,2)-FitResults(1:time,1,2))
xlabel('Time')
ylabel('Position')
text(0,3,'   Peak 1')
text(0,5,'   Peak 2')
text(0,2,'   Separation')
axis([1 100 1 9])