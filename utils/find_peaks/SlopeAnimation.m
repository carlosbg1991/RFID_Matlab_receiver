% Matlab/Octave script visualizing the slope as the tangent to a curve at
% each point. 
x=1:.05:10;
n=length(x);
% y=sin(x); % Creates sine wave
y=exp(-((x-3.5).*2).^2)+.5.*exp(-((x-7.5).*2).^2); % Creates two Gaussian peaks
% y=deriv(y); % Remove comment to see second derivative (slope of first derivative)
plot(x,y)
for j=1:n-1,
    clf
    plot(x,y)
    hold on
    plot(x(j),y(j),'or')
    title('Straight red line is the slope: the tangent to the blue curve at each point')
    xlabel('x')
    ylabel('y')
    P=polyfit([x(j) x(j+1)],[y(j) y(j+1)],1);
    a=P(1);
    b=P(2);
   
    ylim([min(y)-0.1.*min(y) max(y)+0.1.*max(y)])
    plot(x,a.*x+b,'r')
    text(min(x),.1,['   Slope = ' num2str(P(1)) ] );
    drawnow
    pause(.1)
end
