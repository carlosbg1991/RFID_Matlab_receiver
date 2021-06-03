% Demonstration of the reduction of digitization noise by adding random noise
% followed by ensemble averaging. You can change amount of noise in line 6.
format compact
t=0:1000;
S=zeros(1001,100);
Noise=.4;
s=round(10.*sin(.002*pi*t))+Noise.*randn(size(t));
clf;
subplot(2,1,1)
plot(t,s);
axis([0 1000 -12 12])
xlabel('Sine wave with digitization noise casued by truncation/rounding of numbers')
title('Demonstration of ''Digitization Noise''')
subplot(2,1,2)
for k=1:100 % Sets the number of measurement averaged (should be 10 to 100 or so)
  s=round(10.*sin(.002*pi*t)+Noise.*randn(size(t)));
  S(:,k)=s';
end
plot(t,mean(S'))
axis([0 1000 -12 12])
xlabel('Ensemble average of multiple measurements. Try different amounts of noise in line 6')
title(['Effect of adding noise and ensemble averaging. Noise = ' num2str(Noise)] )
e=mean(S')-10.*sin(.002*pi*t);
RMSdifference=sqrt(sum(e.^2));