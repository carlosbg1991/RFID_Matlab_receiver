% Classical Least Squares (CLS) demo (Requires cls.m in the Matlab path)
% Create a noisy model signal of three peaks of known shapes, positions,
% and widths, but unknown heights. Applies Classical Least Squares
% (multilinear regression) using the cls.m function. T. C. O'Haver 2016
format short g
format compact
warning off all
% User settable values
x=200:1:800;
PeakShape=1;
PeakHeights=[5 30 20];
PeakPositions=[400 500 600];
PeakWidths=[70 85 90];
NumPeaks=length(PeakHeights);
extra=0;
Noise=.5;
PeakShift=0; % random peak position and widths shifts between calibration and measurement.
NumTrials=5; % For iterative non-linear peak fit only; higher = slower but more accurate
%
% "start" is the first guess vector for iterative non-linear peak fit
% (length must be 2*NumPeaks). Edit this if you change the number of peaks.
start=[PeakPositions(1) PeakWidths(1) PeakPositions(2) PeakWidths(2) PeakPositions(3) PeakWidths(3)];
%
model=modelpeaks(x,NumPeaks,PeakShape,PeakHeights,PeakPositions,PeakWidths,extra);
y=model+Noise.*randn(size(x)); % Add random noise

% Use the cls.m function to measure the peak heights in the noisy data.
% knowing the NumPeaks,PeakShape,PeakPositions,PeakWidths.
figure(1)
disp(' ')
disp('Figure(1): Classical Least Squares (multilinear regression)')
% Uncomment next two lines to test the effect of random peak position and
% widths shifts between calibration and measurement.
PeakPositions=PeakPositions+PeakShift.*randn(size(PeakPositions));
PeakWidths=PeakWidths+PeakShift.*randn(size(PeakWidths));
tic
MeasuredHeights=cls(x,y,NumPeaks,PeakShape,PeakPositions,PeakWidths,extra);
toc
subplot(2,1,1)
plot(x,y,x,model,'k');
hold on
for peak=1:NumPeaks,
    A(peak,:)= modelpeaks(x,1,PeakShape,PeakHeights(peak),PeakPositions(peak),PeakWidths(peak));
end
plot(x,A')
hold off
title('Classical Least Squares (multilinear regression)')
xlabel('Blue = data       Black = fitted model');ylabel('y')
subplot(2,1,2)
plot(x,y-model,'r');
title('residuals (y-model)')
xlabel('x');ylabel('y')
CLSAccuracy=mean(abs(100.*(PeakHeights-MeasuredHeights)./PeakHeights));
disp([ 'Average peak height accuracy = ' num2str(CLSAccuracy) '%' ] )
disp(' ')

% Use unconstrained iterative non-linear peak fitting using peakfit.m,
% knowing the NumPeaks, PeakShape, and the original unshifted PeakPositions
% and PeakWidths as starting guesses.
figure(2)
disp('Figure(2): Iterative non-linear peak fitting with good starting guesses.')
tic
[Peakfit GOF]=peakfit([x;y],0,0,NumPeaks,PeakShape,extra,NumTrials,start);
toc
PeakfitAccuracy=mean(abs(100.*(PeakHeights-Peakfit(:,3)')./PeakHeights));
disp([ 'Best of ' num2str(NumTrials) ' trial fits.' ] )
disp([ 'Average peak height accuracy = ' num2str(PeakfitAccuracy) '%' ] )
figure(1)