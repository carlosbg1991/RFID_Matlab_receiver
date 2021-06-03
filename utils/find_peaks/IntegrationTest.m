% Exploration of the effect of sampling rate (data density) on the accuracy
% of peak area measurements for single isolated sparsely sampled Gaussian
% peak, y. Computes the percent error in the area calculated by the direct
% sum and trapezoidal method as a function of the number of data points in
% the basewidth (defined as 4 sigma).
figure(1)
clf
clear
trial=1;
for dx=.7:.01:2 % dx is the sampling interval (difference bewtween adjacent x values.
    mu=0;
    sigma=1;
    limits=6; % Number of sigmas plus and minum from peak center mu
    % Define half-Gaussian by defining x for only positive values
    x=mu:dx:limits+mu; % Define sampled points from zero to +limits standard deviations.
    xx=mu-limits:.01:mu+limits; % Make a smooth curve to illustrate Gaussian
    yy=(1/(sqrt(2*pi*sigma^2)).*exp(-((xx-mu).^2)./(2*sigma^2))); % Unit area gaussian. 
    y=(1/(sqrt(2*pi*sigma^2)).*exp(-((x-mu).^2)./(2*sigma^2))); %  Sparsely sampled Gaussian
    ta=trapz(x,y); % Area calculated by trapezoidal method
    sa=sum(y)*dx-((max(y)*dx)/2); % simple half area measurement
    TrueAreaG=.5; % because y was measureed only over half its width
    PercenterrorTrap(trial)=abs(100.*(ta-TrueAreaG)./TrueAreaG); % error of trap method
    PercenterrorS(trial)=abs(100.*(sa-TrueAreaG)./TrueAreaG); % error of simple method
    basewidth=4;
    points(trial)=basewidth/dx; % points = points in basewidth
%    disp('    Trap area   Points/basewidth  % Error T   % Error S   ')
 %   disp([ta,points(trial),PercenterrorTrap(trial),PercenterrorS(trial)] )
    figure(1)
    clf
    
    subplot(2,1,1)
    plot(x,y,'ob',-x,y,'ob',xx,yy,'r')
    xlabel('x')
    ylabel('y')
    title('Red: Gaussian      Blue circles: sampled points ')
    axis([mu-limits limits+mu 0 .5])
    
    subplot(2,1,2)
    plot(points,PercenterrorTrap,'r.',points,PercenterrorS,'b.')
    axis([2 6 0 1])
    title('Effect of data density on the accuracy of  peak area measurements')
    ylabel('Percent error in peak area')
    xlabel('Number of datapoints in basewidth of peak (4*sigma)')
    trial=trial+1;
end
grid