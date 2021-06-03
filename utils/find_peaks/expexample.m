x=1:.2:10; 
a=1;
b=-.9;
y = a.*exp(b.*x);
y=y+y.*.4.*rand(size(x)); 
figure(1)
[coeff,R2]=plotit(x,log(y),1);
ylabel('ln(y)');
title('Plot of x vs the natural log (ln) of y')
aa=exp(coeff(2));
bb=coeff(1);
yy= aa.*exp(bb.*x);
figure(2)
plot(x,y,'r.',x,yy,'g')
xlabel('x');
ylabel('y');
title(['y = a*exp(b*x)     a = ' num2str(aa)  '     b = ' num2str(bb)  '    R2 =  ' num2str(R2) ] ) ;