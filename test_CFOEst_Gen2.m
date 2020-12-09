function fo_est = test_CFOEst_Gen2(varargin)
% test_CFOEst_Gen2 - This script estimates the CFO from an inputted signal
% using Sinewave fitting and modeling. The script loads raw IQ samples from
% an input file and estimates the CFO using a window of samples. The
% boundaries of such window are characterized by the input values 'initFit'
% and 'endFit'. The script can also add synthetic CFO to the inputted sgnal
% using the input parameter 'fo'.
%
% Syntax:  test_CFOEst_Gen2('fileName', 'misc/data/file_source_test_0')
%          test_CFOEst_Gen2('initFit', 4734, 'endFit', 5500)
%
% Inputs:
%    fileName [optional] - string containing the location of the input IQ
%    samples from the RFID reader.
%    initFit [optional] - Beginning of the window of samples used for
%    fitting
%    endFit [optional] - Ending of the window of samples used for
%    fitting
%    fo [optional] - Frequency offset in Hz. This generates a CFO in the
%    input signal of fo Hz.
%
% Outputs: 
%    fo_est - Estimated CFO in Hz
%
%
%------------- BEGIN CODE --------------

%% PARSE INPUTS
p = inputParser;
p.FunctionName = 'test_CFOEst_Gen2';
addParameter(p,'fileName','misc/data/file_source_test', @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'fo',      150, @(x)validateattributes(x,{'numeric'},{'nonempty'}))
addParameter(p,'initFit', 4734, @(x)validateattributes(x,{'numeric'},{'nonempty'}))
addParameter(p,'endFit',  6935, @(x)validateattributes(x,{'numeric'},{'nonempty'}))
parse(p,varargin{:})
fileName = p.Results.fileName;
initFit = p.Results.initFit;
endFit = p.Results.endFit;
fo = p.Results.fo;


addpath('utils/');

%% RX PARAMETERS
[adc_rate, ~, oMF, ~, ~] = fconfig;

%% Input test sequence
% We retrieve the recorded samples without CFO from GNURadio application
fid = fopen(fileName, 'r');
data = fread(fid,inf,'single');
raw_IQ = (data(1:2:end) + 1i.*data(2:2:end)).';

%% CFO
CFO = fo/adc_rate; 
phase_freq_offset = exp(1i*2*pi*(0:numel(raw_IQ)-1)*CFO);
raw_IQ_CFO = raw_IQ .* phase_freq_offset; 

%% Matched Filter
my_filter = ones(1,oMF.n_taps);
output_pulse = conv(my_filter,raw_IQ_CFO);
output_MF = decimate(output_pulse,oMF.decim);

%% Retrieve the portion that we wish to use to estimate the CFO
range_CW = (initFit:endFit);
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
phase_delta = 2*pi*CFO_est_Re*(initFit-1);
s_CFO_corr = offset_Re + (ampl_Re/oMF.n_taps) .* sin(2*pi*(0:N-1)*CFO_est_Re + phase_Re + phase_delta) + ...
             offset_Im + (ampl_Im/oMF.n_taps) .* 1i.*sin(2*pi*(0:N-1)*CFO_est_Im + phase_Im + phase_delta);
         
fo_est = (fo_est_Im+fo_est_Re)/2;

fprintf('Length %d -> CFO_est Real: %.4f (%.4f)\n',N,fo_est_Re,fo);
fprintf('Length %d -> CFO_est Imag: %.4f (%.4f)\n',N,fo_est_Im,fo);