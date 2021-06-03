function y = violetnoise(x)
% Random noise with high-frequency weighted power spectrum with mean zero
% and unit standard deviation. Length equal to the length of x
y=zeros(size(x));
for n=1:2:length(x)-1,
    rn=abs(randn());
    y(n)=rn;
    y(n+1)=-rn;
end