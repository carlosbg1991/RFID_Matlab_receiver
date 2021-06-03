x=1:100000;
% 100000 random numbers using the randn function
randmnums=randn(size(x));
% 100000 pseudo-random numbers using the rand function
fakerandnums=sqrt(3)*(rand(size(x))-rand(size(x))+rand(size(x))-rand(size(x)));
[Nr,Xr]=hist(randmnums);real=peakfit([Xr;Nr]);
[Nf,Xf]=hist(fakerandnums);fake=peakfit([Xf;Nf]);
figure(1);clf;subplot(2,1,1)
realG=real(3).*gaussian([-4:.1:4],real(2),real(4));
plot(Xr,Nr,'o',[-4:.1:4],realG)
title('Normally distributed random numbers using the RANDN function')
xlabel(['Mean = ' num2str(mean(randmnums)) '    Standard deviation = '  num2str(std(randmnums')) ] )
axis([-4           4           0       max(realG)])

figure(1);subplot(2,1,2)
fakeG=fake(3).*gaussian([-4:.1:4],fake(2),fake(4));
plot(Xf,Nf,'o',[-4:.1:4],fakeG)
title('Approximation using using sqrt(3)*(RAND()-RAND()+RAND()-RAND())')
xlabel(['Mean = ' num2str(mean(fakerandnums)) '    Standard deviation = '  num2str(std(fakerandnums')) ] )
axis([-4           4           0       max(fakeG)])
disp('---------------------------------------------------------------')
disp('                  MEAN          STD          Max         Min')
disp(['Real RandNums    ' num2str([mean(randmnums) std(randmnums) max(randmnums)  min(randmnums)]) ] )
disp(['Fake RandNums    ' num2str([mean(fakerandnums) std(fakerandnums) max(fakerandnums)  min(fakerandnums)]) ] )