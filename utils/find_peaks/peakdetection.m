x=1:.1:100;
y=x/100+gaussian(x,30,5)+gaussian(x,70,10);
dy=derivxy(x,y);
plot(x,y,x,dy,'r');
hold on;
plot([30 30],[-.5 max(y)],'c',[70 70],[-.5 max(y)],'c');
hold off
xlabel('X')
ylabel('Y')
title('The locations of the maxima of the peaks (blue) correspond to zero crossings of the first derivative (red) ')
