function [avBER, avPDR] = main_RX_Matlab(varargin)
% GOAL: EMULATE THE RX GNURADIO CODE IN MATLAB

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

%% RX PARAMETERS
[adc_rate, fs, oMF, oG, oTD] = fconfig;

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
N = numel(output_MF);

% %% Retrieve the portion that we wish to use to estimate the CFO
% my_init = 4300;
% my_end = 6098;
% range_CW = (my_init:my_end);
% input_CFOEst = output_MF(range_CW);
% 
% %% Estimate the CFO
% x = (1:numel(range_CW))./fs;
% % REAL
% y = real(input_CFOEst);
% sineParams = CBG_sineFit(x,y);
% offset_Re = sineParams(1);
% ampl_Re = sineParams(2);
% fo_est_Re = sineParams(3);
% phase_Re = sineParams(4);
% % IMAG
% y = imag(input_CFOEst);
% sineParams = CBG_sineFit(x,y);
% offset_Im = sineParams(1);
% ampl_Im = sineParams(2);
% fo_est_Im = sineParams(3);
% phase_Im = sineParams(4);

% %% Correct CFO
% CFO_est_Re = fo_est_Re/fs;
% CFO_est_Im = fo_est_Im/fs;
% phase_delta = 2*pi*CFO_est_Re*(my_init-1);
% s_CFO_corr = offset_Re + (ampl_Re/oMF.n_taps) .* sin(2*pi*(0:N-1)*CFO_est_Re + phase_delta) + ...
%              offset_Im + (ampl_Im/oMF.n_taps) .* 1i.*sin(2*pi*(0:N-1)*CFO_est_Im + phase_delta);
% 
% %% Compare CFO sequences
% s_CFO_decim = decimate([s_CFO zeros(1,oMF.n_taps)],oMF.decim);
% figure;
% subplot(2,1,1); hold on;
% plot((1:numel(s_CFO_decim))./fs.*1e3,real(s_CFO_decim),'linewidth',1.5);
% plot((1:numel(s_CFO_decim))./fs.*1e3,real(s_CFO_corr),'linewidth',1.5);
% ylabel('Real (amplitude)')
% xlabel('Time (ms)');
% legend('original CFO seq','estimated CFO seq');
% set(gca,'FontWeight','bold','fontSize',9);
% subplot(2,1,2); hold on;
% plot((1:numel(s_CFO_decim))./fs.*1e3,imag(s_CFO_decim),'linewidth',1.5);
% plot((1:numel(s_CFO_decim))./fs.*1e3,imag(s_CFO_corr),'linewidth',1.5);
% ylabel('Imaginary (amplitude)')
% xlabel('Time (ms)');
% legend('original CFO seq','estimated CFO seq');
% set(gca,'FontWeight','bold','fontSize',9);        
%          
% output_MF = output_MF.*conj(s_CFO_corr);

%% DISCRETE EVENT SIMULATOR
DES_index = 1;  % place D ES index at the beginning of the file
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