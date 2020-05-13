clear all; 
% close all; 
clc;

addpath('../');

adc_rate = 2e6;  % DO NOT CHANGE
N = 100e4;  % number of samples in the CW

SNR = 30;  % desired SNR in dB
fo = 2e2;  % selected freq offset in Hz

% generate CW
s_CW = ones(1,N);
P_CW = (1/N)*(s_CW*s_CW');

% generate noise
s_n = randn(1,N) + 1i.*randn(1,N);  % noise
P_N = (1/N)*(s_n*s_n');

% apply selected SNR
s_n = sqrt(P_CW/(P_N*db2pow(SNR))) .* s_n;

% add signals
y = s_CW + s_n;

% apply CFO
CFO = fo/adc_rate;
s_CFO = exp(1i*2*pi*(0:N-1)*CFO);
y_CFO = y .* s_CFO;

% apply MF
n_taps = 25;
decim = 5;
my_filter = ones(1,n_taps);
output_pulse = conv(my_filter,y_CFO);
output_MF = decimate(output_pulse,decim);

% estimate CFO using FFT
[fo_est] = f_estimate_CFO(output_MF,adc_rate/decim);

% correct CFO
CFO_est = fo_est/adc_rate;
s_CFO_corr = exp(-1i*2*pi*(0:N-1)*CFO_est);
y_corr = y_CFO.*s_CFO_corr;

figure; hold on;
plot(1:numel(y_CFO),real(y_CFO));
plot(1:numel(y_corr),real(y_corr));
plot(1:numel(s_CFO_corr),real(s_CFO_corr));
ylim([-2 2]);
legend('with CFO','corrected CFO','correction seq');
title(strcat('N=',num2str(N),{' '},'- error of',{' '},num2str(abs(fo_est-fo)),'Hz'))
set(gca,'FontWeight','bold','fontSize',9);

fprintf('N = %d (%d) - Error: %.4f\n',N,round(N/decim),abs(fo_est-fo));