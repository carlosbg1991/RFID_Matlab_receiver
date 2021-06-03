% Demonstration of the linearization of the power transform method by
% taking the nth power of data and then the 1/n power of the measured
% areas. Requires gaussian.m, halfwidth.m and plotit.m downloaded from my
% web site: http://tinyurl.com/cey8rwh

%
% User changable constants
n=1; % power. Must be a positive integer. Power=1 is the "normal" method.
width=3.34; % peak width (FWHM). 
sigma=width/(2*log(2)); % peak width (sigma)
position1=4; % x-axis position of peak 1
position2=8; % x-axis position of peak 2
Resolution=(position2-position1)./(4*sigma); % 4-sigma resolution
noise=.00; % Random noise standard deviation (optional)

% Main processing loop
figure(1)
clf
x=0:.01:12;
for k=1:1:20  
    h(k)=.15*(k/2+2); % Heights of the second peak
    y1=gaussian(x,position1,width); % Peak 1 (fixed height)
    y2=(h(k).*gaussian(x,position2,width));% Peak 2 (variable height)
    a1(k)=trapz(x,y1); % Areas of isolated peak 1
    a2(k)=trapz(x,y2); % Areas of isolated peak 2
    a2p(k)=trapz(x,y2.^n); % Areas of isolated peak 2 raised to power n
    w2p(k)=halfwidth(x,y2.^n);  % Widths of isolated peak 2 raised to power n
    y=y1+y2+noise*randn(size(x)); % Overlapping peak 1 and 2 signal
    yp=y.^n; % Raise to power n
    hold on
    Peak1MaxIndex=val2ind(x,position1); % Index of maximum peak 1
    Peak2MaxIndex=val2ind(x,position2); % Index of maximum peak 2
    valleyy=min(yp(Peak1MaxIndex:Peak2MaxIndex)); % Y value at valley
    valleyx=x(val2ind(yp,valleyy)); % X value at valley
    ValleyIndex=val2ind(x,valleyx); % Index of valley point
    plot(x,yp);
    ap1(k)=sum(yp(1:ValleyIndex))*(x(2)-x(1)); %  perpendicular drop areas of peak 1
    ap2(k)=sum(yp(ValleyIndex:length(x)))*(x(2)-x(1)); %  perpendicular drop areas of peak 2
    apr1(k)=ap1(k).^(1/n); % nth root of perpendicular drop areas of peak 1
    apr2(k)=ap2(k).^(1/n); % nth root of perpendicular drop areas of peak 2
end
text(1,max(yp),[ '4-sigma resolution = ' num2str(Resolution) ] )
text(1,0.95.*max(yp),[ 'Power = ' num2str(n) ] )
text(1,0.9.*max(yp),[ 'Noise = ' num2str(noise) ] )
xlabel('x (time)')
hold off

% Plots response curve
figure(2)
clf
coef=plotit(a2,apr2,1,'or','b-'); xlabel('True area of isolated peak 2, before y is raised to n power') 
ylabel('1/n power of area of signal raised to n power')
title(['Calibration curve of the Power Method. Power = ' num2str(n) '   Resolution = ' num2str(Resolution) ])

