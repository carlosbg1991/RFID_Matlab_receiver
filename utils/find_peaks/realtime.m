% Use of mouse pointer to simulate a signal coming from an instrument or
% sensor. Position your mouse pointer along the y (vertical) axis of the
% graph and click to enter data points as you move the mouse pointer up and
% down the y axis.  Data points are assigned to the vector y (line 17),
% plotted on the graph as black points (line 18), and print out in the
% cvommand window (line 19)
format compact
format short g
clf
maxx=10;
axis([0 10 0 10]);
hold on;
y=zeros(1,maxx);
disp('            n          y')
for n=1:maxx
    [clickX,clickY] = ginput(1);
    y(n)=clickY;
    plot(n,y(n),'k.')
    disp([n y(n)]);
end
hold off
