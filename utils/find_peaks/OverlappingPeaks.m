% Demonstration of gaussfit function for two overlapping peaks
noise=0.001; % Percent random white noise (optional)
x=1:1:100;
% Peak parameters
position1=45;
height1=1;
width1=10;
position2=60;
height2=.2;
width2=10;
% S = signal; N = noise
S=height1.*gaussian(x,position1,width1) + height2.*gaussian(x,position2,width2); % signal vector
N=noise.*randn(size(x)); % noise vector
y=S+N; % Observed signal plus noise

% Delect data ranges for each peak and apply gaussfit function
range1=35:50; % data range for peak 1
range2=57:69; % data range for peak 2
tic
[MeasuredHeight1,MeasuredPosition1,MeasuredWidth1]=gaussfit(x(range1),y(range1));
[MeasuredHeight2,MeasuredPosition2,MeasuredWidth2]=gaussfit(x(range2),y(range2));

figure(1)
clf
plot(x,y,'og')
xlabel('x')
ylabel('y')
title('Applying the gaussfit.m function to overlapping peaks requires carefull selection of data region for each peak')
hold on
plot(x(range1),MeasuredHeight1.*gaussian(x(range1),MeasuredPosition1,MeasuredWidth1),'-r')
plot(x(range2),MeasuredHeight2.*gaussian(x(range2),MeasuredPosition2,MeasuredWidth2),'-b')
hold off
GFtime=toc;
PositionError1=100.*(position1-MeasuredPosition1)./position1;
HeightError1=100.*(height1-MeasuredHeight1)./height1;
WidthError1=100.*(width1-MeasuredWidth1)./width1;
PositionError2=100.*(position2-MeasuredPosition2)./position2;
HeightError2=100.*(height2-MeasuredHeight2)./height2;
WidthError2=100.*(width2-MeasuredWidth2)./width2;
disp(' ')
disp('gaussfit method, applied separately to two peaks. ')
disp('% Parameter errors:   Position   Height     Width')
disp(['Peak 1:              ' num2str(PositionError1) '     '   num2str(HeightError1) '   '  num2str(WidthError1) ])
disp(['Peak 2:              ' num2str(PositionError2) '     '   num2str(HeightError2) '   '  num2str(WidthError2) ])
disp(['Elapsed time, seconds: ' num2str(GFtime)   ])

figure(2)
tic
P=peakfit([x;y],0,0,2,1);
PFtime=toc;
MeasuredPosition1=P(1,2);MeasuredHeight1=P(1,3);MeasuredWidth1=P(1,4);
MeasuredPosition2=P(2,2);MeasuredHeight2=P(2,3);MeasuredWidth2=P(2,4);
PositionError1=100.*(position1-MeasuredPosition1)./position1;
HeightError1=100.*(height1-MeasuredHeight1)./height1;
WidthError1=100.*(width1-MeasuredWidth1)./width1;
PositionError2=100.*(position2-MeasuredPosition2)./position2;
HeightError2=100.*(height2-MeasuredHeight2)./height2;
WidthError2=100.*(width2-MeasuredWidth2)./width2;
disp(' ')
disp('Iterative curve fitting method, using peakfit.m, measures both peaks together')
disp('P=peakfit([x;y],0,0,2,1) ')
disp('% Parameter errors:   Position      Height      Width')
disp(['Peak 1:              ' num2str(PositionError1) '     '   num2str(HeightError1) '   '  num2str(WidthError1) ])
disp(['Peak 2:              ' num2str(PositionError2) '     '   num2str(HeightError2) '   '  num2str(WidthError2) ])
disp(['Elapsed time, seconds: ' num2str(PFtime)   ])
subplot(2,1,1)
title('Iterative curve fitting method, using peakfit.m, measures both peaks together')