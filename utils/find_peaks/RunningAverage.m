% Demonstration of the running average of random numbers with a population
% mean of 1 and a standard deviation of 1.
maxcount=1000;
A=zeros(1,maxcount);
sumy=0;
clf
for cnt=1:maxcount
  y = 1+randn(1);
  sumy=sumy+y;
  A(cnt)=sumy/cnt;
  plot(1:cnt,A(1:cnt),'-')
  axis([1 cnt+1 .5 1.5])
  grid
  xlabel('Count')
  ylabel('Average Value')
  drawnow
end
title('Running average of random numbers with a population mean of 1')