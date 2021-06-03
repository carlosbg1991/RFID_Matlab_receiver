% Demonstration of the linearization of the power transform method by
% taking the 1/n power of the measured areas of data rasied to the power n.
n=5;
figure(1)
width=2;
x=-5:.1:5;
for k=10:-1:1;
    h(k)=k/10;
    y=(h(k).*gaussian(x,0,width)).^n;
    plot(x,y);
    hold on
    a(k)=h(k).*trapz(x,gaussian(x,0,width));
    ap(k)=trapz(x,y).^(1/n);
    w(k)=halfwidth(x,y);
end
hold off

figure(2)
clf
 coef=plotit(a,ap,1);
 xlabel('Peak Area, before being raised to n power')
 ylabel('1/n power of area of signal raised to n power')
 title(['Demonstration of the Power Transform Method. Power = ' num2str(n)])
disp(['Note that the slope of the curve is ' num2str(coef(1)) ' times that of the original power 1 curve.' ])