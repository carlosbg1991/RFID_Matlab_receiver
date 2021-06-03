function mY=medianfilter(y,MedianWidth)
% function mY=medianfilter(y,MedianWidth)
% Performs a median-based filter operation that replaces each elment of y
% with the median of MedianWidth adjacent points (which must be a positive
% integer.  Useful for eliminating narrow spike artifacts in signals.
% Example:
% >> x=-5:.001:5;y=exp(-(x).^2);
% >> for n=1:1000:10000,y(n)=randn(1)+y(n);,end
% >> subplot(1,2,1);plot(x,y);title('Before')
% >> subplot(1,2,2);plot(x,medianfilter(y,1));title('After')
% (c) Tom O'Haver, 2012
mY=y;
MW=abs(round(MedianWidth));
for n=1:length(y)-(1+MW*2),
    mY(n+MW)=median(y((n):(n+1+MW*2)));
end
