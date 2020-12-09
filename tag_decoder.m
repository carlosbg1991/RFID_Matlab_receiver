function [tag_bits,EPC_bits,EPC_hex,EPC_long_hex,h_est,oTD] = tag_decoder(in,oTD)
% [tag_bits,EPC_bits,EPC_hex,oTD] = tag_decoder(in,oTD)
% Performs tasks of fine-grained frame synchronization and FM demodulation
% of RN16 and EPC sequences. For the EPC, it also checks the CRC and
% returns the decoded identifier of the tag in hexadecimal. 
%
% Syntax:  [tag_bits,EPC_bits,EPC_hex,oTD] = tag_decoder(in,oTD)
%
% Inputs:
%    in - IQ samples after the Matched Filter and Gate blocks, representing
%         the coarse synchronization of the packet.
%    oTD - Configuration of the Tag Decoder module.
%
% Outputs:
%    tag_bits - Vector [1 x 16] with the RN16 decoded bits.
%    EPC_bits - Vector [1 x 128] with the decoded EPC bits.
%    EPC_hex - Hexadecimal value of the EPC.
%    h_est - Channel estimation, scalar and complex value.
%    oTD - Updated configuration of the Tag Decoder module.
%
% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
%
% See also: main_RX_Matlab.m

%------------- BEGIN CODE --------------
%% PARAMETERS
tag_bits = zeros(1,oTD.RN16_BITS-1);
EPC_bits = zeros(1,oTD.EPC_BITS-1);
EPC_hex = 65535;  % error code
EPC_long_hex = 65535;  % error code

%% DETECT START OF RN16 SEQUENCE
[start_index, h_est] = tag_sync(in,oTD);

switch oTD.state
    case oTD.SEEK_RN16
        % Store in RN16_samples_complex vector the FM0 modulated RN16_BITS
        number_of_half_bits = 0;
        RN16_samples_complex = [];
        for j = (start_index : oTD.n_samples_TAG_BIT/2 : numel(in))
            number_of_half_bits = number_of_half_bits + 1;
            k = round(j);
            RN16_samples_complex = [RN16_samples_complex in(k)];  %#ok<AGROW>

            if (number_of_half_bits == 2*(oTD.RN16_BITS-1))
                tag_bits = tag_detection_RN16(RN16_samples_complex,h_est,oTD);
                if oTD.PLOT
                    figure(3); plot(real(RN16_samples_complex),imag(RN16_samples_complex),'lineStyle','none','marker','*');
                    xlim([-2 2]); ylim([-2 2]); 
                    hold off;
                end
                break;
            end
        end
        % Update TD state
        oTD.state = oTD.SEEK_EPC;
    case oTD.SEEK_EPC
        % compute RSSI
%         RSSI_dBm = f_compute_RSSI(in, start_index, oTD);
        EPC_bits = tag_detection_EPC(in, start_index, h_est, oTD);
        if(numel(EPC_bits) == oTD.EPC_BITS - 1)
            % float to char -> use Buettner's function
            char_bits = zeros(1,oTD.EPC_BITS-1);
            for i = (1:1:128)
                if (EPC_bits(i)==0)
                    char_bits(i) = 0;
                else
                    char_bits(i) = 1;
                end
            end
            
            result = 0;
            result1 = 0;
            result2 = 0;
            result3 = 0;
            if(check_crc(char_bits,128) == 1)
                for i = (1:1:16)
                    result = result + 2^(16-i)*EPC_bits(96+i);
                end
                for i = (1:1:16)
                    result1 = result1 + 2^(16-i)*EPC_bits(80+i);
                end
                for i = (1:1:16)
                    result2 = result2 + 2^(16-i)*EPC_bits(64+i);
                end
                for i = (1:1:16)
                    result3 = result3 + 2^(16-i)*EPC_bits(48+i);
                end
                EPC_hex = dec2hex(result);
                EPC_long_hex = append( ...
                               append(repmat('0',1,4-length(dec2hex(result3))),dec2hex(result3)), ' ', ...
                               append(repmat('0',1,4-length(dec2hex(result2))),dec2hex(result2)), ' ', ...
                               append(repmat('0',1,4-length(dec2hex(result1))),dec2hex(result1)), ' ', ...
                               append(repmat('0',1,4-length(dec2hex(result))),dec2hex(result)) ...
                               );
            end
        end
        % Update TD state
        oTD.state = oTD.SEEK_RN16;
    otherwise
        error('Wrong TD state');
end
end

function [max_index,h_est] = tag_sync(in,oTD)
    max_index = 0;
    max = 0;

    % Do not have to check entire vector (not optimal)
    for i = (1 : 1 : 1.5*oTD.n_samples_TAG_BIT - 1)
        corr2 = 0;
%         corr = 0;
        % sync after matched filter (equivalent)
        for j = (1 : 1 : 2*oTD.TAG_PREAMBLE_BITS)
            corr2 = corr2 + in( round(i+j*oTD.n_samples_TAG_BIT/2)) * (oTD.TAG_PREAMBLE(j) + 0*1i);
        end
        corr = norm(corr2);
        if (corr > max)
            max = corr;
            max_index = i;
        end
    end

    % Preamble ({1,1,-1,1,-1,-1,1,-1,-1,-1,1,1} 1 2 4 7 11 12))
    h_est = (in(max_index) + in( round(max_index +  1*oTD.n_samples_TAG_BIT/2) ) + ...
                             in( round(max_index +  3*oTD.n_samples_TAG_BIT/2) ) + ...
                             in( round(max_index +  6*oTD.n_samples_TAG_BIT/2) ) + ...
                             in( round(max_index + 10*oTD.n_samples_TAG_BIT/2) ) + ...
                             in( round(max_index + 11*oTD.n_samples_TAG_BIT/2)) ) / ...
                             oTD.TAG_PREAMBLE_BITS;

    % Shifted received waveform by n_samples_TAG_BIT/2
    max_index = max_index + oTD.TAG_PREAMBLE_BITS * oTD.n_samples_TAG_BIT + oTD.n_samples_TAG_BIT/2;
    max_index = max_index + oTD.n_samples_TAG_BIT/2;  % Matlab Correction
end

function tag_bits = tag_detection_RN16(RN16_samples_complex, h_est, oTD)
    % detection + differential decoder (since Tag uses FM0)
    prev = 1;
    tag_bits = zeros(1,oTD.RN16_BITS-1);
    for j = (1 : 1 : numel(RN16_samples_complex)/2)
        result = real( (RN16_samples_complex(2*j-1) - RN16_samples_complex(2*j)) * conj(h_est));

        if (result>0)
            if (prev == 1)
                tag_bits(j) = 0;
            else
                tag_bits(j) = 1;
                prev = 1;
            end
        else
            if (prev == -1)
                tag_bits(j) = 0;
            else
                tag_bits(j) = 1;
                prev = -1;
            end
        end
    end
end

function tag_epc_bits = tag_detection_EPC(in, index, h_est, oTD)
    tag_epc_bits = zeros(1,oTD.EPC_BITS-1);
    magn_squared_samples = abs(in).^2;
    number_steps = 20;
    min_val = oTD.n_samples_TAG_BIT/2.0 - oTD.n_samples_TAG_BIT/2.0/100;
    max_val = oTD.n_samples_TAG_BIT/2.0 + oTD.n_samples_TAG_BIT/2.0/100;

    energy = zeros(1,number_steps);
    for t = (1:1:number_steps)
        for i = (0:1:256)
            energy(t) = energy(t) + ...
                        magn_squared_samples( round(i*(min_val + ...
                           (t-1)*(max_val-min_val)/(number_steps-1)) + index) );
        end
    end
    [~,index_T] = max(energy);
    T =  min_val + index_T*(max_val-min_val)/(number_steps-1);

    % T estimated
%     T_global = T;

    prev = 1;
    temp = zeros(1,2*128);
    for j = (1:1:128)
        temp(2*j-1) = in( round((j-1)*(2*T) + index) );
        temp(2*j)   = in( round((j-1)*(2*T) + T + index) );
         
        result = real((in( round((j-1)*(2*T) + index) ) - ...
                       in( round((j-1)*2*T + T + index) )) * ...
                       conj(h_est) );

        if (result>0)
            if (prev == 1)
                tag_epc_bits(j) = 0;
            else
                tag_epc_bits(j) = 1;
                prev = 1;
            end
        else
            if (prev == -1)
                tag_epc_bits(j) = 0;
            else
                tag_epc_bits(j) = 1;
                prev = -1;
            end
        end
    end
    if oTD.PLOT
        figure(4);
        plot(real(temp),imag(temp),'lineStyle','none','marker','*');
        hold on;
        plot(real(temp.*conj(h_est)),imag(temp.*conj(h_est)),'lineStyle','none','marker','*');
        hold off;
    end
end

function status = check_crc(bits, num_bits)
    num_bytes = num_bits / 8;

    data = zeros(1,num_bytes);
    for i = (1: 1: num_bytes)
        mask = 128;
        data(i) = 0;
        for j = (1: 1: 8)
            if (bits(((i-1) * 8) + j) == 1)
                data(i) = bitor(data(i),mask);
            end
            mask = bitshift(mask,-1);
        end        
    end
    rcvd_crc = bitshift(data(num_bytes - 1),+8) + data(num_bytes);

    crc_16 = 65535;
    for i = (1: 1: num_bytes-2)
        crc_16 = bitxor(crc_16,bitshift(data(i),+8));
        for j = (1: 1: 8)
            if (bitand(crc_16,32768))
                res = dec2bin( bitshift( crc_16, +1 ) );
                if numel(res)>16
                    crc_16 = bin2dec(res(2:end));
                else
                    crc_16 = bin2dec(res);
                end
                crc_16 = bitxor(crc_16,4129);
            else
                res = dec2bin( bitshift( crc_16, +1 ) );
                if numel(res)>16
                    crc_16 = bin2dec(res(2:end));
                else
                    crc_16 = bin2dec(res);
                end
            end
        end
    end
    
    % Complementary
    A = (2^16 - 1) - crc_16;
    bitArray = cell2mat(arrayfun(@dec2bin,A,'un',0));
    crc_16 = bin2dec(bitArray);

    if(rcvd_crc ~= crc_16)
        status = -1;
    else
        status = 1;
    end
end