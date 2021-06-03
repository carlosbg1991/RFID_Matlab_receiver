% Demonstration of the killspikes.m function, which must be in the path.
% function fy=killspikes(x,y,SlopeThreshold,width)
% Related functions: medianfilter.m
% Tom O'Haver 2016 toh@umd.edu
format compact
clf
% Create a smooth peak signal
x=0:.01:10;
y=exp(-(x-5).^2);
% Add spikes of 1,2 or 3 points
for n=1:100:length(x)-100,
    spike=randn(1);
    y(n)=y(n)+spike; % 1 point spikes
    y(n+1)=y(n+1)+spike; % 2 point spikes
    % y(n+2)=y(n+2)+spike; % 3 point spikes
end
y=y+.01.*randn(size(x)); % Add random white noise
% Plot the original spiky data
subplot(2,1,1);plot(x,y);

% Call the killspikes function with selected input arguments
SlopeThreshold=.1;
width=2;
fy=killspikes(x,y,SlopeThreshold,width);
% Compare to median filter
% fy=medianfilter(y,width);

% Plot results
subplot(2,1,2);plot(x,fy)
% figure(2)
% isignal(x,fy);
% figure(3)
% [Fit,Error]=peakfit([x;y]) % Fit the original spiky signal
% figure(4)
% [FitFiltered,ErrorFiltered]=peakfit([x;fy]) % Fit the treated signal