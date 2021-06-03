% testautopeaks.m runs all the examples in the autopeaksplot.m help file,
% additionally plotting the data in figure windows 1 and 2 and numbering
% the peaks. Hint: position the figure windows 1 and 2 side by side for
% viewing. Requires gaussian.m in the path.
format short g
format compact
disp(' ')
disp('Test of autopeaksplot.m')
disp('Example 1:  One input argument; data in single vector (x=1:length(y))')
x=[0:.01:5];y=sin(10*x);
P=autopeaksplot(y);
title('Sine wave, PeakMax = 1')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
pause(1)
disp(' ')

disp('Example 2:  One input argument; data in two columns of a matrix')
x=[0:.01:5]';y=x.*sin(x.^2).^2;M=[x y];
P=autopeaksplot(M);
title('y=x.*sin(x^2)^2')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
pause(1)
disp(' ')

disp('Example 3: Two input arguments; data in separate x and y vectors')
x=[0:.1:100];y=(x.*sin(x)).^2;
P=autopeaksplot(x,y);
title('y=x.*sin(x^2)^2')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
pause(1)

x=[0:.005:2];y=humps(x);
P=autopeaksplot(x,y);
title('Built-in "humps" function')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
pause(1)
disp(' ')

disp('Example 4:  Additional input argument (after the x,y data) to control')
disp('peak sensitivity; higher numbers for more peaks. Average peak height = 10')
x=[0:.1:10];y=5+5.*sin(x)+randn(size(x));
P=autopeaksplot(x,y,3);
title('Example 4: Two broad noisy peaks, average peak height = 10')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
disp(' ')
pause(1)

x=[0:.1:100];y=5+5.*cos(x)+randn(size(x));    
plot(x,y)
P=autopeaksplot(x,y,10);
title('15 narrower noisy peaks, average peak height = 10')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
disp(' ')
pause(1)

x=[0:.1:200];y=5+5.*cos(x)+randn(size(x));
plot(x,y)
P=autopeaksplot(x,y,100);
title('31 even narrower noisy peaks, average peak height = 10')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
pause(1)
disp(' ')

disp('Example 5: Seven input arguments. Specify all peak detections parameters')
disp('Two overlapping noisy peaks, actual peak height = 1')
x=1:.2:100;
y=gaussian(x,40,10)+gaussian(x,50,10)+.01.*randn(size(x));
P=autopeaksplot(x,y,0.00026015,0.031007,19,21,3);
title('Two overlapping noisy peaks, actual peak height = 1')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
pause(1)
disp(' ')

disp('Example 6: Seven input arguments. Specify all peak detections parameters,')
disp('in this case using vectors to optimize for peaks with very different widths.')
disp('Two noisy peaks of different widths, actual peak height = 1')
x=1:.2:100;
y=gaussian(x,20,1.5)+gaussian(x,80,30)+.02.*randn(size(x));
P=autopeaksplot(x,y,[0.001 .0001],[.2 .2],[5 10],[10 100],3);
title('Example 6: Two noisy peaks of different widths, actual peak height = 1')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
pause(1)
disp(' ')

disp('Example 7: Peaks on curved baseline, theoretical height 1. Slight data smoothing. ')
x=-50:.1:50;
y=exp(-(x).^2)+exp(-(x+50*rand()).^2)+.02.*randn(size(x));
y=y+rand().*gaussian(x,-50,100);
y=fastsmooth(y,3,3,1);
P=autopeaksplot(x,y);
title('Example 7: Peaks on curved baseline, theoretical height 1. Slight data smoothing.')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
pause(1)
disp(' ')

y=exp(-(x).^2)+exp(-(x+50*rand()).^2)+.02.*randn(size(x));
y=y+rand().*gaussian(x,-50,100);
y=fastsmooth(y,3,3,1);
P=autopeaksplot(x,y);
title('Example 7: Peaks on curved baseline, theoretical height 1. Slight data smoothing.')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
disp('TEST COMPLETE')