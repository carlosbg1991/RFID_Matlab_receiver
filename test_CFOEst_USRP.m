clear all; 
close all; 
clc;

addpath('../');

N = 1e5;  % number of samples in the CW

% adc_rate = 2e6;  % DO NOT CHANGE
% init = 50;
% v = read_complex_binary('/home/carlos/git/uhd/host/build/examples/usrp_samples.dat',N+init).';

% adc_rate = 200e3;  % DO NOT CHANGE
% init = 1e4;
% v = read_complex_binary ('/home/carlos/git/DATA/B210_noSync_internal.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/B210_sync_octo.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/B210_sync_rfclk_no_amp_1.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/B210_sync_rfclk_with_amp_1.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/B210_s ync_rfclk_with_amp_2.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/B210_sync_rfclk_with_amp_3.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/B210_sync_rfclk_with_amp_4.dat',N+init).';

% new data
adc_rate = 2e6;  % DO NOT CHANGE
init = 10e4;
% v = read_complex_binary('/home/carlos/git/DATA/data_new/B210_internal_no-sync_2Mbps_925Mhz_wireless.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/data_new/B210_octo_2Mbps_925Mhz_wired.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/data_new/B210_octo_2Mbps_925Mhz_wireless.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/data_new/B210_rfclk_2Mbps_925Mhz_wired_1.dat',N+init).';
v = read_complex_binary ('/home/carlos/git/DATA/data_new/B210_rfclk_2Mbps_925Mhz_wireless.dat',N+init).';
% v = read_complex_binary ('/home/carlos/git/DATA/data_new/B210_rfclk_noAmp_2Mbps_900Mhz_wireless_interference.dat',N+init).';

% FFT_rfclk
% adc_rate = 2e6;  % DO NOT CHANGE
% init = 100e4;
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_Amp_900M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_Amp_1800M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_TX2_Amp_860M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_TX2_Amp_900M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_TX2_Amp_940M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_TX2_Amp_980M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_TX2_Amp_flt_860M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_TX2_Amp_flt_900M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX1_TX2_Amp_flt_940M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX2_Amp_940M.dat',N+init).';
% v = read_complex_binary('/home/carlos/git/DATA/FFT_rfclk/usrp_samples_TX2_Amp_1880M.dat',N+init).';

y_CFO = v(init+1:end);  % chop out the first samples

figure; hold on;
plot((1:numel(y_CFO))./adc_rate.*1e3,real(y_CFO),'linewidth',1.5);
plot((1:numel(y_CFO))./adc_rate.*1e3,imag(y_CFO),'linewidth',1.5);
xlabel('Time (ms)');
ylabel('Amplitude');
set(gca,'FontWeight','bold','fontSize',12);

% apply MF
n_taps = 25;
decim = 5;
my_filter = ones(1,n_taps);
output_pulse = conv(my_filter,y_CFO);
output_MF = decimate(output_pulse,decim);

% estimate CFO using FFT
% [fo_est] = f_estimate_CFO(output_MF,adc_rate/decim);
% estimate CFO using Linear Regression
% [estimated_params]=sine_fit(1:numel(output_MF),real(output_MF));
% estimate CFO using Linear Regression v2
% [offset,ampl,phase,f_raw_est] = f_estimate_CFO_LR(1:numel(output_MF),real(output_MF));
% fo_est = f_raw_est/(2*pi);
% estimate CFO using Linear Regression v4
x = (1:numel(output_MF))./(adc_rate/decim);
% REAL
y = real(output_MF);
sineParams = CBG_sineFit(x,y);
offset_Re = sineParams(1);
ampl_Re = sineParams(2);
fo_est_Re = sineParams(3);
phase_Re = sineParams(4);
% IMAG
y = imag(output_MF);
sineParams = CBG_sineFit(x,y);
offset_Im = sineParams(1);
ampl_Im = sineParams(2);
fo_est_Im = sineParams(3);
phase_Im = sineParams(4);

% correct CFO
% CFO_est = fo_est/adc_rate;
% s_CFO_corr = exp(-1i*2*pi*(0:N-1)*CFO_est);
% s_CFO_corr = offset + (ampl/n_taps).*sin(2*pi*(0:N-1)*CFO_est + phase);
% s_CFO_corr = offset + (ampl/n_taps).*exp(1i*2*pi*(0:N-1)*CFO_est + phase);
CFO_est_Re = fo_est_Re/adc_rate;
CFO_est_Im = fo_est_Im/adc_rate;
s_CFO_corr = offset_Re + (ampl_Re/n_taps) .* sin(2*pi*(0:N-1)*CFO_est_Re + phase_Re) + ...
             offset_Im + (ampl_Im/n_taps) .* 1i.*sin(2*pi*(0:N-1)*CFO_est_Im + phase_Im);

%% Compare CFO sequences
figure;
subplot(2,1,1); hold on;
plot((1:N)./adc_rate.*1e3,real(s_CFO_corr),'linewidth',1.5);
ylabel('Real (amplitude)')
xlabel('Time (ms)');
legend('estimated CFO seq');
set(gca,'FontWeight','bold','fontSize',9);
subplot(2,1,2); hold on;
plot((1:N)./adc_rate.*1e3,imag(s_CFO_corr),'linewidth',1.5);
ylabel('Imaginary (amplitude)')
xlabel('Time (ms)');
legend('estimated CFO seq');
set(gca,'FontWeight','bold','fontSize',9);
         
y_corr = y_CFO.*conj(s_CFO_corr);

figure;
subplot(2,1,1); hold on;
plot((1:numel(y_CFO))./adc_rate.*1e3,real(y_CFO),'linewidth',1.5);
plot((1:numel(y_CFO))./adc_rate.*1e3,real(y_corr),'linewidth',1.5);
plot((1:numel(y_CFO))./adc_rate.*1e3,real(s_CFO_corr),'linewidth',1.5);
ylabel('Real (ampl)')
xlabel('Time (ms)');
legend('with CFO','corrected CFO','correction seq');
title(strcat('N=',num2str(N),{' '},', Fs=',num2str(adc_rate.*1e-6),{' '},'MHz, CFO est:',num2str(fo_est_Re),'Hz'));
set(gca,'FontWeight','bold','fontSize',9);
subplot(2,1,2); hold on;
plot((1:numel(y_CFO))./adc_rate.*1e3,imag(y_CFO),'linewidth',1.5);
plot((1:numel(y_CFO))./adc_rate.*1e3,imag(y_corr),'linewidth',1.5);
plot((1:numel(y_CFO))./adc_rate.*1e3,imag(s_CFO_corr),'linewidth',1.5);
ylabel('Imag (ampl)')
xlabel('Time (ms)');
legend('with CFO','corrected CFO','correction seq');
title(strcat('N=',num2str(N),{' '},', Fs=',num2str(adc_rate.*1e-6),{' '},'MHz, CFO est:',num2str(fo_est_Im),'Hz'));
set(gca,'FontWeight','bold','fontSize',9);