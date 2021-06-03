clf
x=0:.003:10;
y=sin(210.*x)+sin(250.*x);
for n=1:2:3300,
  if randn()>2,
    spike=rand();
    y(n)=spike+y(n);
  end
end
y=y+.01.*randn(size(y));
% isignal([x;y],4.998,9.999,3,10,0,0,0,10,1000,0,0);
y1=fastsmooth(y,9,3);
subplot(2,1,1);
plot(x,y);
title('Original Signal')
subplot(2,1,2);
plot(x,y1)
title('Spikes extracted')
% To detect peaks, use iPeak:
% P=ipeak([x;y1],0,0.01,2e-005,1,5,5,2,0);
% To observe the frequency spectrum of the original signal:
% isignal([x;y],4.998,9.999,0,3,0,0,0,10,1000,0,0,1);