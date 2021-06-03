clf
x=0:.003:10;
y=sin(-(x).^2);
for n=1:3300,
  if randn()>2,
    y(n)=rand()+y(n);
  end
end
y=y+.01.*randn(size(y));
y1=abs(fastsmooth((deriv2(y)).^2,3,2));
subplot(2,1,1);
plot(x,y);
title('Original Signal')
subplot(2,1,2);
plot(x,y1)
title('Spikes extracted')
% To detect peaks, use iPeak:
% P=ipeak([x;y1],0,0.1,2e-005,1,3,3,0.2,0);
% To observe the frequency spectrum of the original signal:
% isignal([x;y],5,10,0,3,0,0,0,0,0,0,0,1);