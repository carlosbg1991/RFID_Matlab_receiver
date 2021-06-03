% TimeTrial is a script that compares the execution 
% times for eighteen signal processing computations 
% running on different versions of Matlab and on Octave,
% using the Matlab/Octave script TimeTrial.m. Must have
% all referenced scripts and functions in the path.
%
% Hardware: standard desktop PC (Intel Core 
% i5, 3 Ghz) running Windows 10 home. Matlab 
% Mobile: iPad, IOS 12.1, with Internet 
% connection via 5G Wi-Fi.
% 
% The main conclusions that can be drawn 
% from all these results are:
% 
% (1) All the versions of Matlab, including 
% the online and iPad ones, provide
% computational speeds that are mostly with 
% a factor of two of each other, and 
% 
% (2) Octave is much slower than any Matlab 
% version.

disp('Conventional "for loop" vs vector/matrix notation (Matlab timing) Loop:')
k=1000;B=zeros(k,k);tic;x=1:k;for n=1:k;B(:,n)=n.*x;end;toc1=toc
disp('Matlab 2009 elapsed time is 1.446651 seconds.')
disp('Matlab 2017b: .006 seconds.')
disp('Matlab Online, R2018b Update 1: 0.003 seconds.')

disp('1')
disp('Matlab/Octave vector/matrix notation:')
k=1000;clear A;tic;x=1:k;A=[1:k]'*x;toc2=toc
disp('Matlab 2009 elapsed time is 0.004306 seconds.')
disp('Matlab 2017b: 0.0015 seconds.')
disp('Matlab Online, R2018b Update 1:  0.0068 seconds.')
disp('time improvement factor for vector/matrix notation:')
disp(toc1/toc2)

disp('2')
disp('Random numbers')
tic;RANDtoRANDN;t1=toc
disp('Matlab 2009: 0.11 seconds.')
disp('Matlab 2017b: 4.7 seconds.')
disp('Matlab Online, R2018b Update 1: 0.36 seconds.')
disp('Octave 2.1 seconds.')

disp('3')
disp('Comparison of types of noise') 
tic;noisetest;toc
disp('Matlab 2009: 6.5 seconds.')
disp('Matlab 2017b: 4.3 seconds.')
disp('Matlab Online, R2018b) Update 1): 3.996692 seconds.')
disp('Octave 35 seconds.')

disp('4')
disp('Smoothing and plotting')
tic;plot(fastsmooth(randn(1,100000),1000,3,0));toc
disp('Matlab 2009: 0.015 seconds.')
disp('Matlab 2017b Home: 0.013 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.011 seconds.')
disp('Octave 3 seconds.')

disp('5')
disp('Median filter')
tic;x=-5:.01:5;y=exp(-(x).^2);
for n=1:100:length(x),y(n)=randn(1)+y(n);end
subplot(1,2,1);plot(x,y);title('Before');
subplot(1,2,2);plot(x,medianfilter(y,1));title('After');toc
disp('Matlab 2009: 0.04 seconds.')
disp('Matlab 2017b Home: 0.043 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.056 seconds.')
disp('Octave 12 seconds.')

disp('6')
disp('smoothed second derivative using ProcessSignal')
tic;x=[0:.01:10]';y=10./(1+(x/4).^2)+exp(-(x-5).^2)+.1.*randn(size(x));
subplot(2,1,1);plot(x,y);
subplot(2,1,2);plot(x,ProcessSignal(x,y,2,60,3,0,0,0,0,0,0));toc
disp('Matlab 2009: .03 seconds.')
disp('Matlab 2017b Home: 0.03 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.023 seconds.')
disp('Octave .18 seconds.')

disp('7')
disp('Quantitative analysis by derivative spectroscopy')
tic;DerivativeDemo;toc
disp('Matlab 2009: 0.06 seconds.')
disp('Matlab 2017b Home: 0.2 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.23 seconds.')
disp('Octave 3 seconds.')

disp('8')
disp('Fourier Filter')
tic;plot(FouFilter(tan(1:100000),1,19,2,2,0));toc
disp('Matlab 2009: 0.02 seconds.')
disp('Matlab 2017b Home: 0.012 seconds.')
disp('Matlab Online, R2018b) Update 1):  seconds.')
disp('Octave 0.35 seconds.')

disp('9')
disp('Polynomial least squares and plotting')
tic;x=0:100;y=100+(x*100)+100.*randn(size(x));[coef,RSquared,~]=plotit(x,y,1);toc
disp('Matlab 2009: 0.0067 seconds.')
disp('Matlab 2017b Home: 0.03 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.054 seconds.')
disp('Octave 0.984 seconds.')

disp('10')
disp('Quadratic Least squares fit to Gaussian peak')
tic;x=50:.1:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
[Height,Position,Width]=gaussfit(x,y);
plot(x,y,'o',x,Height.*gaussian(x,Position,Width));toc
disp('Matlab 2009: 0.0024 seconds.')
disp('Matlab 2017b Home: 0.023 seconds.')
disp('Matlab Online, R2018b) Update 1):  0.024 seconds.')
disp('Octave 0.3 seconds.')

disp('11')
disp('Multilinear Regression')
tic;RegressionDemo2;toc
disp('Matlab 2009: 0.015423 seconds.')
disp('Matlab 2017b Home: 0.079676 seconds.')
disp('Matlab Online, R2018b) Update 1):  0.1567 seconds.')
disp('Octave 0.29689 seconds.')

disp('12 ')
disp('findpeaks function on noisy synthetic data')
tic;DemoFindPeak;toc
disp('Matlab 2009: 0.12 seconds.')
disp('Matlab 2017b Home: .04 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.060 seconds.')
disp('Octave 1.5 seconds.')

disp('13 ')
disp('Demo Find Peak SNR')
tic;DemoFindPeakSNR;toc
disp('Matlab 2009: 0.16 seconds.')
disp('Matlab 2017b Home: .05 seconds.')
disp('Matlab Online, R2018b) Update 1): .045 seconds.')
disp('Octave 2.5 seconds.')

disp('14 ')
disp('Peakfit.m with plotting enabled')
tic;x=0:.1:10;y=exp(-(x-5).^2);peakfit([x' y']);toc
disp('Matlab 2009: 0.03 seconds.')
disp('Matlab 2017b Home: 0.06 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.11 seconds.')
disp('Octave 0.44 seconds.')

disp('15 ')
disp('Peakfit.m with plotting and bootstrap error estimates')
tic;x=0:.1:10;y=10*exp(-(x-4).^2)+randn(size(x));
[FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit([x' y'],0,0,1);toc
disp('Matlab 2009: 0.7 seconds.')
disp('Matlab 2017b Home: 1.1 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.91 seconds.')
disp('Octave 1.8 seconds.')

disp('16 ')
disp('first-order least-squares fit error propagation')
tic;LinearFiMC;toc
disp('Matlab 2009: 0.574 seconds.')
disp('Matlab 2017b Home: 0.4 seconds.')
disp('Matlab Online, R2018b) Update 1): 0.9 seconds.')
disp('Octave 1.9 seconds.')

disp('17 ')
disp('first-order least-squares fit error propagation by bootstrap')
tic;TestLinearFit;toc
disp('Matlab 2009: 0.44 seconds.')
disp('Matlab 2017b Home: 3.4 seconds.')
disp('Matlab Online, R2018b) Update 1): 3.1 seconds.')
disp('Octave 1.9 seconds.')

disp('18 ')
disp('Monte Carlo comparison of full-peak iterative fit with Gaussfit method.')
tic;GaussFitMC2(100,100,100,50,1,1000);toc
disp('Matlab 2009: 7 seconds.')
disp('Matlab 2017b Home: 5.6 seconds.')
disp('Matlab Online, R2018b) Update 1): 4.7 seconds.')
disp('Octave 11 seconds.')