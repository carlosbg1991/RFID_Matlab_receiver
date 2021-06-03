format compact
format short g
Deltat=1; % Sampling interval (1/sample rate)
t=[0:14]'; % time axis
frequency=.1; % Frequency
y=sin(2*pi*frequency*t); % Create a sine wave of that frequency
clf
subplot(2,1,1); % Plot and label it in upper half of figure window
plot(t,y);
title('Signal')
xlabel('time, t')

fy=fft(y); % Compute Fourier transform
sy=fy .* conj(fy); % Multiply by complex conjugate
plotrange=1:length(fy)/2; % Define plotrange as half length of y
ps=real(sy(plotrange)); % ps = power spectrum
f=(plotrange-1)./range(t); % Create frequency axis (reciprocal of t)
subplot(2,1,2); % Plot and label it in lower half of figure window
plot(f,ps)
title('Power spectrum')
xlabel('frequency, f')
[Height,Position,Width]=lorentzfit(f,ps')