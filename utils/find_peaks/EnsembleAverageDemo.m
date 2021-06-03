% Demonstration of ensemble averaging with very low signal-to-noise ratios.
% Change MaxPeakHeight (line 8) to make the problem harder or easier.
increment=0.1;% Change increment to change the number of points in signal.
x=-100:increment:300; 
Noise=2; % Change if you wish.
Baseline=0; % Change if you wish.
MaxPeakHeight=1; % Change if you wish. SNR range will be from zero to this value.
PeakPosition=100; % Change if you wish.
PeakWidth=100;% Change if you wish.
SignalToNoiseRatio=MaxPeakHeight./Noise;
figure(1)
clf
AverageY=zeros(size(x));
% Change number of signals to average in the line below
for trial=1:500,
    % Compute simulated noisy signal
    % You can change 'gaussian' to 'lorentzian' or 'rectanglepulse' or
    % any other shape in the line below
    y=Baseline+MaxPeakHeight.*gaussian(x,PeakPosition,PeakWidth)+Noise.*randn(size(x));
    % Compute ensemble average
    AverageY=AverageY+y;
    % Plot original signal and average
    subplot(1,2,1);plot(x,y,'.');
    xlabel('x')
    ylabel('y')
    title('Original signal')
    subplot(1,2,2);plot(x,AverageY./trial,'.')
    xlabel('x')
    ylabel('total y')
    title('Ensemble averaged signal')
    drawnow
    % pause(.1) % to slow it down
end