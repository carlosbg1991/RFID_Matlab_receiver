% Matlab/Octave script OnePeakOrTwo.m creates a signal that might be
% interpreted as either one peak at x=3 on a curved baseline or as two
% peaks at x=.5 and x=3, depending on context. In this demo, the
% findpeaksG.m function used called twice, with two different values of
% SlopeThreshold to demonstrate.
x=0:.01:5;
y=2.*gaussian(x,.5,5)+gaussian(x,3,1)+.02.*randn(size(x));
figure(1)
clf
plot(x,y,'b.');
disp('One peak at x=3 on a curved baseline, ignoring the baseline peak at x=0.5')
disp('findpeaksG(x,y,0.0001,0.1,30,30,3)')
disp('    Peak#    Position    Height    Width    Area')
disp(findpeaksG(x,y,0.0001,0.1,30,30,3))
xlabel('X')
ylabel('Y')
title('Does this signal have two peaks or one peak on a curved baseline?')

disp('Same signal, two peaks to be detected, at x=0.5 and x-3.')
disp('Reduce SlopeThreshold to .00001')
disp('findpeaksG(x,y,0.00001,0.1,30,30,3)')
disp('    Peak#    Position    Height    Width    Area')
disp(findpeaksG(x,y,0.00001,0.1,30,30,3))

