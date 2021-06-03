format short g

% Generate pulse signal
t = 0:10:10000; % Sampling frequency
duty=.5; % Percent duty cycle (50 = square wave)
y = 1400+1300.*(sign(duty+sin(2*pi*.001*t))+.02.*randn(size(t)));
plot(t,y);
xlabel('Seconds');
ylabel('Amplitude');

S=findsquarepulse(t,y,1000);
disp('     Pulse number  Start time    Height        Width ')
disp(S)