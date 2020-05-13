% EPC_bits = [0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,...
%             0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,1,...
%             1,0,1,1,0,0,1,0,1,1,0,1,1,1,0,1,1,1,0,1,...
%             1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,...
%             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,...
%             0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,1,...
%             0,0,1,1,1,1,1,0];
EPC_bits = [0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,...
            0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,...
            0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,...
            0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,...
            1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,0,1,1,1,...
            0,0,0,0,1,1,0,1,1,0,1,0,1,0,0,1,1,1,1,1,...
            0,1,1,0,0,1,0,1];

status = check_crc(EPC_bits, numel(EPC_bits));
         
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