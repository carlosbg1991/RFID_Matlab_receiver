% Plotting separate graphs using 'subplot' and 'title'

x=0:.01:6;
clear A

A=[1:10]'*x;
subplot(2,2,1);
plot(x,A)
title('subplot(2,2,1) showing y=n*x')
xlabel('x')
ylabel('y')

A=[1:10]'*sin(x);
subplot(2,2,2);
plot(x,A)
title('subplot(2,2,2) showing y=n*sin(x)')
xlabel('x')
ylabel('y')

clear A
A=[1:10]'*sqrt(x);
subplot(2,2,3);
plot(x,A)
title('subplot(2,2,3) showing y=n*sqrt(x)')
xlabel('x')
ylabel('y')

clear A
A=[1:10]'*x.^2;
subplot(2,2,4);
plot(x,A)
title('subplot(2,2,4) showing y=n*x^2')
xlabel('x')
ylabel('y')