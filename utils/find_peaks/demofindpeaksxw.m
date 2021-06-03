% Find, plot, and number noisy peaks with unknown positions
x=-10:.2:50;
y=exp(-(x).^2)+exp(-(x-50*rand()).^2)+exp(-(x-50*rand()).^2);
y=y+.02.*randn(size(x));
plot(x,y,'m.')
P=findpeaksxw(x,y,0.001,0.2,1,2,3);
text(P(:,2),P(:,3),num2str(P(:,1)))
disp('           peak #    Position    Height       Width')
disp(P)