% Demonstration of the running average of random numbers with a population
% mean of 1 and a standard deviation of 1.
clear
maxcount=30;
S=zeros(1,maxcount);
sumy=0;
clf
for cnt=2:maxcount
    y = randn(1,maxcount);
    S(cnt)=std(y(1:cnt));
end
plot(1:cnt,S,'-')
grid
xlabel('Count')
ylabel('Standard deviation')
title('Running standard deviation of random numbers with a population standard deviation of 1')