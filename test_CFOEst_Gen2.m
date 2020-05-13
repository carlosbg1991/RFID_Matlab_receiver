clear all;
close all;
clc;

addpath('utils/');

%% RX PARAMETERS
[adc_rate, fs, oMF, oG, oTD] = fconfig;

%% Input test sequence
% We retrieve the recorded samples without CFO from GNURadio application
fid = fopen('../../../misc/data/file_source_test', 'r');
data = fread(fid,inf,'single');
raw_IQ = (data(1:2:end) + 1i.*data(2:2:end)).';

%% CFO
fo = 150;
CFO = fo/adc_rate; 
phase_freq_offset = exp(1i*2*pi*(0:numel(raw_IQ)-1)*CFO);
raw_IQ_CFO = raw_IQ .* phase_freq_offset; 

%% Matched Filter
my_filter = ones(1,oMF.n_taps);
output_pulse = conv(my_filter,raw_IQ_CFO);
output_MF = decimate(output_pulse,oMF.decim);

%% Retrieve the portion that we wish to use to estimate the CFO
my_init = 4734;
my_end = 6935;
range_CW = (my_init:my_end);
N = numel(range_CW);
input_CFOEst = output_MF(range_CW);

%% Estimate the CFO
x = (1:numel(input_CFOEst))./(adc_rate/oMF.decim);
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

%% Correct CFO
CFO_est_Re = fo_est_Re/adc_rate;
CFO_est_Im = fo_est_Im/adc_rate;
phase_delta = 2*pi*CFO_est_Re*(my_init-1);
s_CFO_corr = offset_Re + (ampl_Re/oMF.n_taps) .* sin(2*pi*(0:N-1)*CFO_est_Re + phase_Re + phase_delta) + ...
             offset_Im + (ampl_Im/oMF.n_taps) .* 1i.*sin(2*pi*(0:N-1)*CFO_est_Im + phase_Im + phase_delta);

fprintf('Length %d -> CFO_est Real: %.4f (%.4f)\n',N,fo_est_Re,fo);
fprintf('Length %d -> CFO_est Imag: %.4f (%.4f)\n',N,fo_est_Im,fo);