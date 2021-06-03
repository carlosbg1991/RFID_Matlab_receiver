clear
format compact
start=1;  % Starting value in first year
rate=.07; % Annual rate of return
Noise=.5;
NumTrials=30;
disp(['TrueRateOfReturn = ' num2str(rate)])
x=1:66; % time in years
y=start.*(1+rate).^x;
for trial=1:NumTrials,
    % Add proportional random noise to y data
    yy=y+Noise.*y.*randn(size(x));
    %  Method 1: fit noisy log y transfomed data to a straight line
    figure(1)
    [coeff,trR2]=plotit(x,log(abs(yy)),1);
    start1(trial)=exp(coeff(2));
    rate1(trial)=exp(coeff(1))-1;
    xlabel('Time, years')
    ylabel('Natural log of value')
    title('Straight-line fit to log(value) transformed data')
    % Method 2: use iterative line-linear fit to noisy x,y data
    figure(2)
    [rate2(trial),start2(trial),FittingError,fsR2]=fitshape1(x,yy,0);
    xlabel('Time, years')
    ylabel('Value')
    title('Iterative fit to y=start*(1+rate)^x to noisy x,y data')
end
meanstart1=mean(start1);stdstart1=std(start1);rsdstart1=stdstart1./meanstart1;
meanstart2=mean(start2);stdstart2=std(start2);rsdstart2=stdstart2./meanstart2;
meanrate1=mean(rate1);stdrate1=std(rate1);rsdrate1=100.*(stdrate1./meanrate1);
meanrate2=mean(rate2);stdrate2=std(rate2);rsdrate2=100.*(stdrate2./meanrate2);
range1=max(rate1)-min(rate1)
range2=max(rate2)-min(rate2)
text(min(x),max(y)-.05.*(max(y)-min(y)),['   Rate of return = ' num2str(meanrate2)] );
disp('                           Rate          % RSD    high-low range')
disp(['Coordinate transformation: ' num2str([meanrate1 rsdrate1 range1])])
disp(['Iterative curve fitting:   ' num2str([meanrate2 rsdrate2 range2])])
disp(' ')