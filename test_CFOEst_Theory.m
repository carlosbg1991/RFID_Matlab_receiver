 function [CFO_error] = test_CFOEst_Theory(varargin)

%% Check correct number of input arguments. If none, then set-up default
% parameters from configuration functions
if (nargin==4)
    if ~ (isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
            isnumeric(varargin{3}) && isnumeric(varargin{4}))
        error('ERROR: Input parameters should be numeric\n');
    end
    
    N = varargin{1};
    SNR = varargin{2};
    fo = varargin{3};
    po = varargin{4};
elseif (nargin==0)
    clear; close all;
    % Pre-configured initial values
    N = 1e5;  % number of samples in the CW
    SNR = 30;  % desired SNR in dB
    fo = 102;  % selected freq offset in Hz
    po = 2*pi*rand(1);  % phase offset
end

addpath('utils/');  % include utils functions

adc_rate = 2e6;  % DO NOT CHANGE

%% generate CWcd 
s_CW = ones(1,N); 
P_CW = (1/N)*(s_CW*s_CW');

%% generate noise
s_n = randn(1,N) + 1i.*randn(1,N);  % noise
P_N = (1/N)*(s_n*s_n');

%% apply selected SNR 0
s_n = sqrt(P_CW/(P_N*db2pow(SNR))) .* s_n;

%% add signals
y = s_CW + s_n;

%% apply CFO
CFO = fo/adc_rate;  % cycles per sample
s_CFO = exp(1i*(2*pi*(0:N-1)*CFO + po));
y_CFO = y .* s_CFO;

%% apply MF
n_taps = 25;
decim = 5;
my_filter = ones(1,n_taps);
output_pulse = conv(my_filter,y_CFO);
output_MF = decimate(output_pulse,decim);

%% estimate CFO using Linear Regression
% Retrieve the portion used to compute the CFO
% my_init = 0.5*(1/CFO);
my_init = 1;
my_end = numel(output_MF);
input_CFOEst = output_MF(my_init:my_end);
x = (1:numel(input_CFOEst))./(adc_rate/decim);
% REAL
y = real(input_CFOEst);
sineParams = CBG_sineFit(x,y);
offset_Re = sineParams(1);
ampl_Re = sineParams(2);
fo_est_Re = sineParams(3);
phase_Re = sineParams(4);
% IMAG
y = imag(input_CFOEst);
sineParams = CBG_sineFit(x,y);
offset_Im = sineParams(1);
ampl_Im = sineParams(2);
fo_est_Im = sineParams(3);
phase_Im = sineParams(4);

%% CFO correction
CFO_est_Re = fo_est_Re/adc_rate;
CFO_est_Im = fo_est_Im/adc_rate;
phase_delta = 2*pi*CFO_est_Re*(my_init-1);
s_CFO_corr = offset_Re + (ampl_Re/n_taps) .* sin(2*pi*(0:N-1)*CFO_est_Re + phase_Re + phase_delta) + ...
             offset_Im + (ampl_Im/n_taps) .* 1i.*sin(2*pi*(0:N-1)*CFO_est_Im + phase_Im + phase_delta);

%% CFO error
CFO_error = min(abs(fo_est_Re-fo),abs(fo_est_Im-fo));
         
%% Compare CFO sequences
% figure;
% subplot(2,1,1); hold on;
% plot((1:numel(s_CFO))./adc_rate.*1e3,real(s_CFO),'linewidth',1.5);
% plot((1:numel(s_CFO))./adc_rate.*1e3,real(s_CFO_corr),'linewidth',1.5);
% ylabel('Real (amplitude)')
% xlabel('Time (ms)');
% legend('original CFO seq','estimated CFO seq');
% set(gca,'FontWeight','bold','fontSize',9);
% subplot(2,1,2); hold on;
% plot((1:numel(s_CFO))./adc_rate.*1e3,imag(s_CFO),'linewidth',1.5);
% plot((1:numel(s_CFO))./adc_rate.*1e3,imag(s_CFO_corr),'linewidth',1.5);
% ylabel('Imaginary (amplitude)')
% xlabel('Time (ms)');
% legend('original CFO seq','estimated CFO seq');
% set(gca,'FontWeight','bold','fontSize',9);
% 
% y_corr = y_CFO.*conj(s_CFO_corr);
% 
% figure;
% subplot(2,1,1); hold on;
% plot((1:numel(y_CFO))./adc_rate.*1e3,real(y_CFO),'linewidth',1.5);
% plot((1:numel(y_CFO))./adc_rate.*1e3,real(y_corr),'linewidth',1.5);
% plot((1:numel(y_CFO))./adc_rate.*1e3,real(s_CFO_corr),'linewidth',1.5);
% ylabel('Real (ampl)')
% xlabel('Time (ms)');
% legend('with CFO','corrected CFO','correction seq');
% title(strcat('N=',num2str(N),{' '},', Fs=',num2str(adc_rate),{' '},' error=',num2str(abs(fo-fo_est_Re)),'Hz'));
% set(gca,'FontWeight','bold','fontSize',9);
% subplot(2,1,2); hold on;
% plot((1:numel(y_CFO))./adc_rate.*1e3,imag(y_CFO),'linewidth',1.5);
% plot((1:numel(y_CFO))./adc_rate.*1e3,imag(y_corr),'linewidth',1.5);
% plot((1:numel(y_CFO))./adc_rate.*1e3,imag(s_CFO_corr),'linewidth',1.5);
% ylabel('Imag (ampl)')
% xlabel('Time (ms)');
% legend('with CFO','corrected CFO','correction seq');
% title(strcat('N=',num2str(N),{' '},', Fs=',num2str(adc_rate),{' '},' error=',num2str(abs(fo-fo_est_Im)),'Hz'));
% set(gca,'FontWeight','bold','fontSize',9);