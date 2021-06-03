function fy=killspikes(x,y,Threshold,width)
% function fy=killspikes(x,y,Threshold,width)
% Function for removing positive and negative spikes from time series data.
% Each time it finds a positive or negative jump in the data between x(n)
% and x(n+1) that exceeds "Threshold", it replaces the next "width" points
% of data with a linearly interpolated segment spanning x(n) to x(n+width).
% Example: 
% x=1:9;y=[1 2 3 4 10 6 7 8 9];
% subplot(2,1,1);plot(x,y,'.');
% fy=killspikes(x,y,1,2);
% subplot(2,1,2);plot(x,fy,'.')
% See killspikesdemo.m
% Related function: medianfilter.m
% Tom O'Haver 2016 toh@umd.edu
fy=y;
for j=2:length(y)-3,
    if abs(y(j)-y(j+1)) > Threshold, % if difference between points is larger than SlopeThreshold
        for k=1:width,
            fy(j+k-1)=interp1([x(j-1) x(j+1+width)],[fy(j-1) fy(j+1+width)],x(j+k-1));
        end
    end
end
