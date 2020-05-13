function [CFO_est] = f_estimate_CFO(input_samples,Fs)
%% FFT            
T = 1/Fs;               % Sampling period       
L = numel(input_samples);  % Length of signal
t = (0:L-1)*T;          % Time vector

t_new = (0:1:L-1)*T;      % Time vector
input_samples = interp1(t,input_samples,t_new,'linear');

Y = fft(input_samples);

P2 = abs(Y/L);  % normalize FFT
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

[~,idx] = max(P1);
CFO_est = f(idx);

figure; hold on;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

end