function BER = f_compute_BER(EPC_bits, EPC_bits_ground_truth)
    BER = sum(abs(EPC_bits_ground_truth-EPC_bits))./numel(EPC_bits);
end