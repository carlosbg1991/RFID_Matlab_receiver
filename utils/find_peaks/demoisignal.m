% Self-running demo script of iSignal 6
% Must have isignal.m, peakfit.m, and plotit.m in the path.
pausetime=2; % Determines speed of demo
x=0:.005:2;y=humps(x);Data=[x;y];
disp('The built-in Matlab Humps function, displayed in iSignal')
isignal(Data,0.6,4,0,0,1,0,0,0,0,0,0);
subplot(2,1,1)
title('The built-in Matlab Humps function, displayed in iSignal')
pause(pausetime)
isignal(Data,0.245,0.615,4,3,1,0,0,0,0,0,0,0);
subplot(2,1,1)
title('Viewing the first peak near x=0.4')
drawnow
pause(pausetime)
isignal(Data,0.945,0.735,4,3,1,0,0,0,0,0,0);
subplot(2,1,1)
title('Viewing the second peak near x=0.9')
pause(pausetime)
isignal(Data,0.9,0.5,1,3,1,4);
subplot(2,1,1)
title('4th derivative of the peak at x=0.9')
pause(pausetime)
isignal(Data,0.3,0.5,1,3,1,0,1,220,5400);
subplot(2,1,1)
title('Peak sharpening applied to the peak at x=0.3')
pause(pausetime)
clf
x=-5:.01:5;
y=exp(-(x).^2);
for n=1:1000,
    if randn()>2,y(n)=rand()+y(n);
    end
end
isignal([x;y],0,20,0,3,0,0,0,10,1000,0,0,0);
subplot(2,1,1)
title('Smooth signal with many sharp spikes')
pause(pausetime)
isignal([x;y],0,20,0,3,0,0,0,10,1000,0,1,0);
subplot(2,1,1)
title('Spikes removed using median filter (M key)')
pause(pausetime)
x=0:.1:60; y=sin(x)+sin(10.*x);
[pY,PowerSpectrum]=isignal([x;y],30,30,4,3,1,0,0,1,0,0,0,1);
pause(pausetime)
clf
plot(PowerSpectrum(:,1),PowerSpectrum(:,2))
grid
title('Power spectrum plotted on linear coordinates')
pause(pausetime)
loglog(PowerSpectrum(:,1),PowerSpectrum(:,2))
grid
title('Power spectrum plotted on log-log coordinates')
pause(pausetime)
% Demonstration script for iSignal function. It generates a test signal
% consisting of 4 peaks on a curved baseline, then runs iSignal
%   T. C. O'Haver, December 2013, January 2015
increment=.1;
x=1:increment:700;
% For each simulated peak, compute the amplitude, position, and width
% Peak 1 is the baseline
amp=[100 .5 2 3 3.5];  % Amplitudes of the peaks  (CHANGE if desired)
pos=[-200 100 250 400 600];   % Positions of the peaks (CHANGE if desired)
wid=[700 50 50 50 50];   % Widths of the peaks (CHANGE if desired)
Noise=.02; % Amount of random noise added to the signal. (CHANGE if desired) 
% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      A(k,:)=gaussian(x,pos(k),wid(k)); % Gaussian or Lorentzian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
DataMatrix=[x' y']; % Assembles x and y vectors into data matrix
disp('-----------------------------------------------------------------')
disp('Detection and measurement of 4 peaks heights 1, 2, 3, 3.5, with ')
disp('equal widths, superimposed in a very strong curved baseline.')
disp('The objective is to extract a measure that is proportional to ')
disp('the peak height but independent of the baseline strength. In this demo, ')
disp('differentiation (with smoothing) is used to supress the baseline.')
isignal(DataMatrix);
disp('Signal consisting of 4 weak peaks on a strong curved baseline.')
subplot(211);
title('Signal consisting of 4 weak peaks on a strong curved baseline.')
pause(pausetime)
isignal(DataMatrix,250.5,499,0,1,0,0,0,10,1000,0,0,0);
disp('Press Ctrl-A to zoom out to see entire signal')
subplot(211);
title('Press Ctrl-A to zoom out to see entire signal')
pause(pausetime)
isignal(DataMatrix,250.5,499,0,1,0,1,0,10,1000,0,0,0);
disp('Press D key to compute first derivative')
subplot(211);
title('Press D key to compute first derivative')
pause(pausetime)
isignal(DataMatrix,250.5,499,0,1,0,2,0,10,1000,0,0,0);
disp('Press D key again to compute second derivative')
subplot(211);
title('Press D key again to compute second derivative')
pause(pausetime)
isignal(DataMatrix,250.5,499,1,12,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 12 point rectangular smooth')
pause(pausetime)
isignal(DataMatrix,250.5,499,1,25,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 25 point rectangular smooth')
pause(pausetime)
isignal(DataMatrix,250.5,499,1,51,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 51 point rectangular smooth')
pause(pausetime)
isignal(DataMatrix,250.5,499,1,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point rectangular smooth')
pause(pausetime)
isignal(DataMatrix,250.5,499,2,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point triangular smooth')
pause(pausetime)
isignal(DataMatrix,250.5,499,3,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point Gaussian smooth')
pause(pausetime)
isignal(DataMatrix,250.5,499,4,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point Savitsky-Golay smooth')
pause(pausetime)
isignal(DataMatrix,250.5,499,3,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point Gaussian smooth')
pause(pausetime)
[pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(DataMatrix,100,147.4,3,95,0,2,0,10,1000,0,0,0);
s1=maxy-miny;
subplot(211);
title('Measure peak-to-peak signal of peak 1')
pause(pausetime)
[pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(DataMatrix,250,147.4,3,95,0,2,0,10,1000,0,0,0);
s2=maxy-miny;
subplot(211);
title('Measure peak-to-peak signal of peak 2')
pause(pausetime)
[pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(DataMatrix,400,147.4,3,95,0,2,0,10,1000,0,0,0);
s3=maxy-miny;
subplot(211);
title('Measure peak-to-peak signal of peak 3')
pause(pausetime)
[pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(DataMatrix,600,147.4,3,95,0,2,0,10,1000,0,0,0);
s4=maxy-miny;
subplot(211);
title('Measure peak-to-peak signal of peak 4')
pause(pausetime)
clf
plotit([amp(2) amp(3) amp(4) amp(5)],[s1 s2 s3 s4],1,'or');    
title('Calibration plot of derivative signal ampliude against true peak amplitude')  
xlabel('True amplitude of the 4 original peaks')
ylabel('Peak-to-peak derivative amplitude of each peak')

% % Optional: Compare to peakfit method using 2-peak fits, one peak for
% % the background and the other for each peak at x=100, 250, 400, 600)
% % Measurement of peak heights
% figure(2)
% peakshape=1; % 1=Gaussian. Also try different "incorrect" peakshape numbers
% [FitResults,FitError]=peakfit(DataMatrix,97,190,2,peakshape,0,10, [99.86066      49.63798     -201.1655      701.4819]);
% s1=min(FitResults(:,3));
% drawnow
% [FitResults,FitError]=peakfit(DataMatrix,239.7,190,2,peakshape,0,10, [249.925      50.17847     -208.6503      706.8251]);
% s2=min(FitResults(:,3));
% drawnow
% [FitResults,FitError]=peakfit(DataMatrix,392,190,2,peakshape,0,10, [399.9098      50.59616     -286.9422      748.8739]);
% s3=min(FitResults(:,3));
% drawnow
% [FitResults,FitError]=peakfit(DataMatrix,591.8,190,2,peakshape,0,10, [599.9833      49.98055     -210.0941      704.7717]);
% s4=min(FitResults(:,3));
% drawnow
% figure(3)
% clf
% plotit([amp(2) amp(3) amp(4) amp(5)],[s1 s2 s3 s4],1,'or');  
% title('Accuracy of peak height measurements')
% 
% % Measurement of peak positions
% figure(4)
% peakshape=1; % 1=Gaussian. Also try different "incorrect" peakshape numbers
% [FitResults,FitError]=peakfit(DataMatrix,97,190,2,peakshape,0,10, [99.86066      49.63798     -201.1655      701.4819]);
% p1=max(FitResults(:,2));
% drawnow
% [FitResults,FitError]=peakfit(DataMatrix,239.7,190,2,peakshape,0,10, [249.925      50.17847     -208.6503      706.8251]);
% p2=max(FitResults(:,2));
% drawnow
% [FitResults,FitError]=peakfit(DataMatrix,392,190,2,peakshape,0,10, [399.9098      50.59616     -286.9422      748.8739]);
% p3=max(FitResults(:,2));
% drawnow
% [FitResults,FitError]=peakfit(DataMatrix,591.8,190,2,peakshape,0,10, [599.9833      49.98055     -210.0941      704.7717]);
% p4=max(FitResults(:,2));
% drawnow
% figure(5)
% clf
% plotit([pos(2) pos(3) pos(4) pos(5)],[p1 p2 p3 p4],1,'or'); 
% title('Accuracy of peak posiiton measurements')


disp('End of demo.')