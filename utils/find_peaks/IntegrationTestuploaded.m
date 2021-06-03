% Exploration of the effect of sampling rate (data density) on the accuracy
% of peak area measurements for single isolated sparsely sampled Gaussian
% peak, y. Computes the percent error in the area calculated by the direct
% sum and trapezoidal method as a function of the number of data points in
% the basewidth (defined as 4 sigma).
format compact
format short g
figure(1)
clf
clear
trial=1;
disp('      Points      T area    % Error T      S area     % Error S   ')
for dx=1:.01:2 % dx is the sampling interval (difference bewtween adjacent x values.
    mu=0;
    sigma=1;
    startx=-5;
    x=mu-5:dx:5+mu; % Define sampled points from -5 to +5 standard deviations.
    xx=startx+mu:.01:startx+10+mu; % Make a smooth curve to illustrate Gaussian
    yy=(1/(sqrt(2*pi*sigma^2)).*exp(-((xx-mu).^2)./(2*sigma^2))); % Unit area gaussian. 
    y=(1/(sqrt(2*pi*sigma^2)).*exp(-((x-mu).^2)./(2*sigma^2))); %  Sparsely sampled Gaussian
    ta=trapz(x,y); % Area calculated by trapezoidal method
    sa=sum(y)*dx; % simple area measurement
    TrueAreaG=1; % 
    PercenterrorTrap(trial)=abs(100.*(ta-TrueAreaG)./TrueAreaG); % error of trap method
    PercenterrorS(trial)=abs(100.*(sa-TrueAreaG)./TrueAreaG); % error of simple method
    basewidth=4;
    points(trial)=basewidth/dx; % points = points in basewidth
    disp([points(trial),ta,PercenterrorTrap(trial),sa, PercenterrorS(trial)] )
    figure(1)
    clf
    
    subplot(2,1,1)
    plot(x,y,'ob',xx,yy)
    xlabel('x')
    ylabel('y')
    title('Red: Gaussian      Blue circles: sampled points ')
    axis([startx+mu startx+10+mu 0 .5])
    
    subplot(2,1,2)
    plot(points,PercenterrorTrap,'r.',points,PercenterrorS,'b.')
    axis([2 4 0 1])
    title('Effect of data density on the accuracy of  peak area measurements')
    ylabel('Percent error in peak area')
    xlabel('Number of datapoints in basewidth of peak (4*sigma)')
    trial=trial+1;
end

