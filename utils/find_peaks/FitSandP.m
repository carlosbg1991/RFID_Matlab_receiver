% load SandPfrom1950

range1=1:8304; 
range2=8305:16608;
T=data(:,8)';
V=data(:,5)';

% Method 1: use iterative line-linear fit to x,y data
figure(1)
[rate1,start1,FittingError,fsR2]=fitshape1(T,V,0);
xlabel('Time, years')
ylabel('Value')
title('Iterative fit to y=start*(1+rate)^x to noisy x,y data')
text(min(T),max(V)-.05.*(max(V)-min(V)),['   Rate of return = ' num2str(rate1)] );

% Method 2: fit log y transfomed data to a straight line
figure(2)
[coeff,trR2]=plotit(T,log(V),1);
start2=exp(coeff(2));
rate2=exp(coeff(1))-1;
xlabel('Time, years')
ylabel('Natural log of value')
title('Straight-line fit to log(value) transformed data')

disp(['Iterative curve fitting:   ' num2str(rate1)])
disp(['Coordinate transformation: ' num2str(rate2)])

