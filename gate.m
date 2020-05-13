function [state,start_index,oG] = gate(in,oG)
% [state,start_index,oG] = gate(in,oG)
% Performs tasks of coarse-grained frame synchronization by looking at the
% pulses percieved from the Query and ACK packets sent from Reader.
%
% Syntax:  [tag_bits,EPC_bits,EPC_hex,oTD] = tag_decoder(in,oTD)
%
% Inputs:
%    in - IQ samples after the Matched Filter and Gate blocks, representing
%         the coarse synchronization of the packet.
%    oG - Configuration of the Gate module.
%
% Outputs:
%    state - 1 if Query/ACK is found, 0 otherwise.
%    start_index - Start index of the RN16 and EPc packets representing the
%                  coarse estimation of the start of those packets.
%    oG - Updated configuration of the Gate module.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% 
% See also: main_RX_Matlab.m

%------------- BEGIN CODE --------------
%% PARAMETERS
N = numel(in);
NUM_PULSES_COMMAND = 5;
THRESH_FRACTION = 0.75; 
NEG_EDGE = 0;
POS_EDGE = 1;

signal_state = NEG_EDGE;
state = 0;  % Initialize state to NOT_FOUND
start_index = 0;  % Points to the start of the RN16

% oG.n_samples = 0;
% oG.num_pulses = 0;

for i = 1:N
    %% TRACK AMPLITUDE
    sample_ampl = abs(in(i));
    oG.avg_ampl = oG.avg_ampl + (sample_ampl - oG.win_samples(oG.win_index)) / oG.win_length;
    oG.win_samples(oG.win_index) = sample_ampl;
    oG.win_index = mod((oG.win_index),oG.win_length) + 1;
    
    % Threshold for detecting negative/positive edges
    sample_thresh = oG.avg_ampl * THRESH_FRACTION;

    %% TRACK DC OFFSET (only during T1)
    oG.dc_est =  oG.dc_est + (in(i) - oG.dc_samples(oG.dc_index)) / oG.dc_length;
    oG.dc_samples(oG.dc_index) = in(i);
    oG.dc_index = mod((oG.dc_index),oG.dc_length) + 1;
    
    oG.n_samples = oG.n_samples + 1;
    
    % Potitive edge -> Negative edge
    if( sample_ampl < sample_thresh && signal_state == POS_EDGE)
        oG.n_samples = 0;
        signal_state = NEG_EDGE;
    % Negative edge -> Positive edge
    elseif (sample_ampl > sample_thresh && signal_state == NEG_EDGE)
        signal_state = POS_EDGE;
        if (oG.n_samples > oG.n_samples_PW/2)
            oG.num_pulses = oG.num_pulses + 1;
        else
            oG.num_pulses = 0;
            oG.n_samples = 0;
        end
    end

    if(oG.n_samples==oG.n_samples_T1 && signal_state==POS_EDGE && oG.num_pulses>NUM_PULSES_COMMAND)
        start_index = i;
        state = 1;
        
        oG.n_samples = 0;
        oG.num_pulses = 0;
        
        return;
    end
end


end

%% EOF