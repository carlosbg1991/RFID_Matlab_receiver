% testautopeaks.m runs all the examples in the autopeaks.m help
% file, additionally plotting the data and numbering the peaks
format short g
format compact
disp(' ')
disp('Test of autopeaks.m')
disp('Example 1:  One input argument; data in single vector (x=1:length(y))')
x=[0:.01:5];y=sin(10*x);
P=autopeaks(y);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
plot(1:length(y),y)
title('Sine wave, PeakMax = 1')
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 2:  One input argument; data in two columns of a matrix')
x=[0:.01:5]';y=x.*sin(x.^2).^2;M=[x y];
P=autopeaks(M);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
plot(x,y)
title('y=x.*sin(x^2)^2')
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 3: Two input arguments; data in separate x and y vectors')
x=[0:.1:100];y=(x.*sin(x)).^2;
P=autopeaks(x,y);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
plot(x,y)
title('y=x.*sin(x^2)^2')
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

x=[0:.005:2];y=humps(x);
P=autopeaks(x,y);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
plot(x,y)
title('Built-in "humps" function')
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 4:  Additional input argument (after the x,y data) to control')
disp('peak sensitivity; higher numbers for more peaks:')
x=[0:.1:10];y=5+5.*sin(x)+randn(size(x));
plot(x,y)
y=fastsmooth(y,3,3,1);
P=autopeaks(x,y,3);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
title('Example 4: Two broad noisy peaks, average peak height = 10')
text(P(:,2),P(:,3),num2str(P(:,1)))
disp(' ')
pause(1)

x=[0:.1:100];y=5+5.*cos(x)+randn(size(x));    
plot(x,y)
P=autopeaks(x,y,10);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
plot(x,y)
title('15 narrower noisy peaks, average peak height = 10. ')
text(P(:,2),P(:,3),num2str(P(:,1)))
disp(' ')
pause(1)

x=[0:.1:200];y=5+5.*cos(x)+randn(size(x));
plot(x,y)
title('Even narrower noisy peaks, average peak height = 10')
P=autopeaks(x,y,100);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 5: Seven input arguments. Specify all peak detections parameters')
x=1:.2:100;
y=gaussian(x,40,10)+gaussian(x,50,10)+.01.*randn(size(x));
plot(x,y,'c.')
title('Example 5: Two overlapping noisy peaks, actual peak height = 1')
P=autopeaks(x,y,0.00026015,0.031007,19,21,3);
text(P(:,2),P(:,3),num2str(P(:,1)))
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
pause(1)
disp(' ')

disp('Example 6: Seven input arguments. Specify all peak detections parameters,')
disp('in this case using vectors to optimize for peaks with very different widths.')
x=1:.2:100;
y=gaussian(x,20,1.5)+gaussian(x,80,30)+.02.*randn(size(x));
P=autopeaks(x,y,[0.001 .0001],[.2 .2],[5 10],[10 100],3);
plot(x,y,'c.')
title('Example 6: Two noisy peaks of different widths, actual peak height = 1')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(P)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 7: Find, measure, and plot noisy peaks with unknown positions')
x=-50:.2:50;
y=exp(-(x).^2)+exp(-(x+50*rand()).^2)+.02.*randn(size(x));
y=y+rand().*gaussian(x,-50,100);
plot(x,y,'m')
title('Example 7: Peaks on curved baseline, theoretical height 1. ')
P=autopeaks(x,y);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
text(P(:,2),P(:,3),num2str(P(:,1)))
disp(P)
pause(1)
disp(' ')
y=exp(-(x).^2)+exp(-(x+50*rand()).^2)+.02.*randn(size(x));
y=y+rand().*gaussian(x,-50,100);
plot(x,y,'m')
title('Example 7: Peaks on curved baseline, theoretical height 1. ')
P=autopeaks(x,y);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
text(P(:,2),P(:,3),num2str(P(:,1)))
disp(P)
disp('TEST COMPLETE')