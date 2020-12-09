function [adc_rate, fs, oMF, oG, oTD] = fconfig(varargin)

%% PARSE INPUTS
p = inputParser;
validPlot = @(x) islogical(x);
defaultPlot = false;
addOptional(p,'PLOT',defaultPlot,validPlot);
parse(p,varargin{:});
PLOT = p.Results.PLOT;

%% GLOBAL CONFIG
adc_rate = 2e6;  % sampling rate at the RFIC of the SDR

%% MATCHED FILTER (MF)
oMF.n_taps = 25;
oMF.decim = 5;  % decimation factor for the Matched filter
fs = round(adc_rate / oMF.decim);  % sampling rate in the digital system

%% GATE (G) control params
oG.G_window_size = 800;
oG.win_length = 250 * fs / 1e6;
oG.dc_length  = 120 * fs / 1e6;
oG.n_samples_T1 = 240 * fs / 1e6;
oG.n_samples_PW = 12 * fs / 1e6;
oG.avg_ampl = 0;
oG.dc_est = 0;
oG.win_index = 1;
oG.dc_index = 1;
oG.win_samples = zeros(1,oG.win_length);
oG.dc_samples = zeros(1,oG.dc_length);
oG.n_samples = 0;
oG.num_pulses = 0;

%% TAG_DECODER (TD) control params
oTD.T_READER_FREQ = 40e3;
oTD.TAG_BIT_D   = 1 / oTD.T_READER_FREQ * 1e6;
oTD.n_samples_TAG_BIT = oTD.TAG_BIT_D * fs / 1e6;
oTD.TAG_PREAMBLE_BITS = 6;
oTD.EPC_BITS = 129;
oTD.TAG_PREAMBLE = [1,1,0,1,0,0,1,0,0,0,1,1];
oTD.RN16_BITS = 17;
oTD.SEEK_RN16 = 1;
oTD.SEEK_EPC = 2;
oTD.state = oTD.SEEK_RN16;
oTD.to_copy(1) = (oTD.RN16_BITS + oTD.TAG_PREAMBLE_BITS)*oTD.n_samples_TAG_BIT + ...
                   2*oTD.n_samples_TAG_BIT;  % for RN16
oTD.to_copy(2) = (oTD.EPC_BITS + oTD.TAG_PREAMBLE_BITS)*oTD.n_samples_TAG_BIT + ...
                  2*oTD.n_samples_TAG_BIT + oTD.n_samples_TAG_BIT;  % for EPC

oTD.PLOT = PLOT;

              
              
              
% EOF