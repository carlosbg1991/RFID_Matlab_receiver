% testautofindpeaks.m runs all the examples in the autofindpeaks.m help
% file, additionally plotting the data and numbering the peaks (like
% autofindpeaksplot.m)
format short g
format compact
disp(' ')
disp('Test of autofindpeaks.m')
disp('Example 1:  One input argument; data in single vector (x=1:length(y))')
x=[0:.01:5];y=sin(10*x);
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(y);
plot(1:length(y),y)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 2:  One input argument; data in two columns of a matrix')
x=[0:.01:5]';y=x.*sin(x.^2).^2;M=[x y];
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(M);
plot(x,y)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 3: Two input arguments; data in separate x and y vectors')
x=[0:.1:100];y=(x.*sin(x)).^2;
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(x,y);
plot(x,y)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
x=[0:.005:2];y=humps(x);
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(x,y);
disp(P)
plot(x,y)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 4:  Additional input argument (after the x,y data) to control')
disp('peak sensitivity; higher numbers for more peaks:')
x=[0:.1:10];y=5+5.*sin(x)+randn(size(x));
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(x,y,3);
plot(x,y)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
x=[0:.1:100];y=5+5.*cos(x)+randn(size(x));    
plot(x,y)
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(x,y,10);
plot(x,y)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
x=[0:.1:1000];y=5+5.*cos(x)+randn(size(x));
plot(x,y)
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(x,y,100);
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 5: Seven input arguments. Specify all peak detections parameters')
x=1:.2:100;
y=gaussian(x,40,10)+gaussian(x,50,10)+.01.*randn(size(x));
plot(x,y,'c.')
P=autofindpeaks(x,y,0.00026015,0.031007,19,21,3);
text(P(:,2),P(:,3),num2str(P(:,1)))
disp(P)
pause(1)
disp(' ')

disp('Example 6: Seven input arguments. Specify all peak detections parameters,')
disp('in this case using vectors to optimize for peaks with very different widths.')
x=1:.2:100;
y=gaussian(x,20,1.5)+gaussian(x,80,30)+.02.*randn(size(x));
plot(x,y,'c.')
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(x,y,[0.001 .0001],[.2 .2],[5 10],[10 100],3);
disp(P)
text(P(:,2),P(:,3),num2str(P(:,1)))
pause(1)
disp(' ')

disp('Example 7: Find, measure, and plot noisy peaks with unknown positions')
x=-50:.2:50;
y=exp(-(x).^2)+exp(-(x+50*rand()).^2)+.02.*randn(size(x));
plot(x,y,'m')
disp('           peak #    Position      Height     Width        Area')
P=autofindpeaks(x,y);
text(P(:,2),P(:,3),num2str(P(:,1)))
disp('           peak #    Position    Height       Width        Area')
disp(P)
