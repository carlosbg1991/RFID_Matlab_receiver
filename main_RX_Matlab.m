function [avBER, avPDR] = main_RX_Matlab(varargin)
% main_RX_Matlab - This script emulates a RFID Receiver. It takes the raw
% input samples and decodes them to ultimately print the unique identifier
% by terminal. 
%
% This code extends from the following GNURadio code:
% https://github.com/nkargas/Gen2-UHF-RFID-Reader
% To generate new traces, just run the reader and store them in a binary
% file. Then, call this script with those recorded samples. Note that this
% is a passive receiver, so it does not interact with the tags, it merely
% reads whatever it has been transmitted by the reader and tags.
%
% The also outputs BER and EVM measurements based on prior knowledge of the
% tag in range. For instance, one could run the RFID GNURadio reader with
% only one tag in range and store its EPC (in bits) under the variable
% EPC_bits_ground_truth.
%
% In addition, this code allows to set the frequency offset as though
% the transmitter and receiver antennas were missing frequency 
% synchronization.
%
% Syntax:  main_RX_Matlab('fileName','misc/data/file_source_test_0')
%          main_RX_Matlab('fo',250)  % frequency offset in Hz
%          main_RX_Matlab('PLOT',250)  % plot IQ samples
%          main_RX_Matlab('LOGS',true)  % activate decoding LOGS
%
% Inputs:
%    fileName [optional] - string containing the location of the input IQ
%    samples from the RFID reader.
%    fo [optional] - Frequency offset in Hz.
%    PLOT [optional] - Flag (True or False) to plot the intermediate steps
%    in the decoder. This flag is intended for debugging purposes as it
%    plots the detected frames and also the decoded symbols in each slot.
%    LOGS [optional] - Flag (True or False) to print the logs
%    sf [optional] - Spurious additive signal (frequency in Hz).
%    sa [optional] - Spurious additive signal (amplitude in linear scale).
%
% Outputs:
%    avBER - Average measured BER
%    avPDR - Average measured PER
%
%
%------------- BEGIN CODE --------------
persistent p;

if isempty(p)
    p = inputParser;
    p.FunctionName = 'main_RX_Matlab';
    addParameter(p,'fileName','misc/data/file_source_test_0', ...
        @(x)validateattributes(x,{'char'},{'nonempty'}))
    addParameter(p,'fo',0, ...
        @(x)validateattributes(x,{'numeric'},{'nonempty'}))
    addParameter(p,'PLOT',false, ...
        @(x)validateattributes(x,{'logical'},{'nonempty'}))
    addParameter(p,'LOGS',false, ...
        @(x)validateattributes(x,{'logical'},{'nonempty'}))
    addParameter(p,'sf',2.2e3, ...
        @(x)validateattributes(x,{'numeric'},{'nonempty'}))
    addParameter(p,'sa',0.6 * 0.0 / 25, ...
        @(x)validateattributes(x,{'numeric'},{'nonempty'}))
end

parse(p,varargin{:})

sa = p.Results.sa;
fileName = p.Results.fileName;
PLOT = p.Results.PLOT;
LOGS = p.Results.LOGS;
fo = p.Results.fo;
sf = p.Results.sf;

%% CONFIGURE PATH
addpath('utils/');

%% RX PARAMETERS
[adc_rate, fs, oMF, oG, oTD] = fconfig(PLOT);

%% Input test sequence
fid = fopen(fileName, 'r');
data = fread(fid,inf,'single'); 
raw_IQ = (data(1:2:end) + 1i.*data(2:2:end)).';
% EPC_bits_ground_truth = [0  0  1  1  0  0  0  0  0  0  0  0  0  0  ...
%                          0  0  0  0  1  1  0  0  0  0  0  0  0  0  ...
%                          1  0  0  0  0  0  1  1  0  0  1  1  1  0  ...
%                          1  1  0  0  1  0  1  1  0  1  1  1  0  1  ...
%                          1  1  0  1  1  0  0  1  0  0  0  0  0  0  ...
%                          0  1  0  1  0  0  0  0  0  0  0  0  0  0  ...
%                          0  0  0  0  0  0  0  0  0  0  0  0  0  0  ...
%                          0  0  0  0  0  0  0  0  1  0  0  1  1  1  ...
%                          0  1  1  0  1  1  0  1  0  0  1  1  1  1  1  0];
EPC_bits_ground_truth = [0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,...
                         0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,...
                         0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,...
                         0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,1,...
                         0,0,1,0,0,1,1,1,0,0,0,0,1,1,0,1,1,0,1,0,1,0,0,...
                         1,1,1,1,1,0,1,1,0,0,1,0,1];

%% CFO
CFO = fo/adc_rate;
s_CFO = exp(1i*2*pi*(0:numel(raw_IQ)-1)*CFO);
raw_IQ_CFO = raw_IQ .* s_CFO;

%% Spurious additive signal
SAF = sf/adc_rate;
s_SAS = sa.*exp(1i*2*pi*(0:numel(raw_IQ)-1)*SAF);
raw_IQ_CFO = raw_IQ_CFO + s_SAS;

%% Matched Filter
my_filter = ones(1,oMF.n_taps);
output_pulse = conv(my_filter,raw_IQ_CFO);
output_MF = decimate(output_pulse,oMF.decim);

% output = f_CFOEst_alt(output_MF, oMF, fs, CFO);

%% DISCRETE EVENT SIMULATOR
DES_index = 1;  % place DES index at the beginning of the file
DES_index_max = length(output_MF);  % DES cannot go over this index
num_decoded = 0;
avBER = 0;
avPDR = 1;
while(DES_index_max >= DES_index + 2*max(oTD.to_copy))
    %% GATE 
    % Pre-process samples for GATE
    input_G = output_MF(DES_index:DES_index+oG.G_window_size);
    if PLOT
        figure(1); plot(1:numel(input_G),real(input_G),1:numel(input_G),imag(input_G));
    end
    [state_G,start_idx,oG] = gate(input_G,oG);
    
    if state_G == 1
        %% TAG DECODER
        % Pre-process samples for TAG DECODER
        DES_index = DES_index + start_idx;
        % Remove DC component
        TD_window_size = oTD.to_copy(oTD.state);
        input_TD = output_MF(DES_index:DES_index+TD_window_size) - ...
                   oG.dc_est;
        if PLOT
            figure(2); plot(1:numel(input_TD),real(input_TD),1:numel(input_TD),imag(input_TD));
        end
        [RN16_bits,EPC_bits,EPC_hex,EPC_long_hex,h_est,oTD] = tag_decoder(input_TD,oTD);
        
        if oTD.state == oTD.SEEK_EPC
            if LOGS
                fprintf('%d -> %s\n',num_decoded,num2str(RN16_bits));
            end
        else
            BER = f_compute_BER(EPC_bits, EPC_bits_ground_truth);
            PDR = f_compute_PDR(BER);
            avBER = (num_decoded/(num_decoded+1))*avBER + ...
                    1/(num_decoded+1)*BER;
            avPDR = (num_decoded/(num_decoded+1))*avPDR + ...
                    1/(num_decoded+1)*PDR;
            if LOGS
                fprintf('%d -> H_est: %.3f (%.3f)\n',num_decoded,abs(h_est),angle(h_est));
                fprintf('%d -> %s\n',num_decoded,num2str(EPC_hex));
                fprintf('%d -> %s\n',num_decoded,EPC_long_hex);
                fprintf('%d -> BER of %.4f\n',num_decoded,BER);
                fprintf('%d -> EV: %s\n',num_decoded,num2str(abs(EPC_bits - EPC_bits_ground_truth)));
            end
            num_decoded = num_decoded + 1;
        end
        
        % Update 
        DES_index = DES_index + TD_window_size;
    else
        DES_index = DES_index + oG.G_window_size;
    end
end

%% PRINT RESULTS
% fprintf ('%.4f -> BER: %.4f  PDR: %.4f\n', sa, avBER, avPDR);
fprintf ('%.4f %.4f\n', sa, avPDR);



%EOF