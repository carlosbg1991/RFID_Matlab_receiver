function STDError=SmallSampleSTD(SampleSize)
% Demonstraton of the Law of Large Numbers. Computes the sample standard
% deviation of small samples of numbers (SampleSize) drawn from a
% population of normally distributed random numbers with a true population
% standard deviation of 1.00. Plots the histogram of 1,000,000 repeat
% samples of size 'SampleSize'. Returns the relative percent deviation from
% the true population standard deviation.
figure(1)
NumRepeats=1000000;
y = randn(SampleSize,NumRepeats);
STDy=std(y);
[H,edges]=histcounts(STDy);
plot(edges(1:length(H)),H);
PeakSTD=edges(val2ind(H,max(H)));
xlabel('Distribution of the standard deviations of samples')
title(['Sample size = ' num2str(SampleSize)   ])
STDError=100.*(PeakSTD-1);

function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is
% closest to val If more than one element is equally close, returns vectors
% of indicies and values Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);