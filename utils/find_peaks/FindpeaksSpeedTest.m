% Testing template for findpeaks functions
x=[0:.1:1000]';
y=x.*(1+cos(x)).^2+randn(size(x));
clf
plot(x,y)
tic;
[PKS,LOCS]= findpeaks(y,'MINPEAKHEIGHT',1,'MINPEAKDISTANCE',31);
elapsedtimeSPT=toc;
NumPeaksSPT=length(PKS);
PeaksPerSecondSPT=round(NumPeaksSPT./elapsedtimeSPT);

tic;
Px=findpeaksx(x,y,.0001,1,25,2,3);
elapsedtimex=toc;
NumPeaksx=max(Px(:,1));
PeaksPerSecondx=round(NumPeaksx./elapsedtimex);

tic;
Pg=findpeaksG(x,y,.0001,1,50,60,2);
elapsedtimeg=toc;
NumPeaksg=max(Pg(:,1));
PeaksPerSecondg=round(NumPeaksg./elapsedtimeg);
% 
% tic;
% Pi=ipeak([x,y],0,1,3.9992e-006,14,23,12,20,0);
% elapsedtimeiPeak=toc;
% NumPeaksiPeak=max(Pi(:,1))
% PeaksPerSecond=NumPeaksiPeak./elapsedtimeiPeak
disp('--------------------------------------------------------------------')
disp('Function         Number of peaks    Elapsed time    Peaks per second')
disp(['findpeaks (SPT)       ' num2str(NumPeaksSPT) '            ' num2str(elapsedtimeSPT) '           ' num2str(PeaksPerSecondSPT)  ])
disp(['findpeaksx            ' num2str(NumPeaksx) '            ' num2str(elapsedtimex) '         ' num2str(PeaksPerSecondx)  ])
disp(['findpeaksG            ' num2str(NumPeaksg) '            ' num2str(elapsedtimeg) '          ' num2str(PeaksPerSecondg)  ])