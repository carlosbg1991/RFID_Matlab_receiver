function power_dbm = compute_RSSI(in, ini_bits, oTD)
threshold = 0.85;
Zc = 50;  % value for impedance of USRP
average = 0;

% Samples to compute over (RN16 vs EPC)
if (oTD.state == oTD.SEEK_RN16)
    n_samples_bits = oTD.RN16_BITS * oTD.n_samples_TAG_BIT;
else
    n_samples_bits = oTD.EPC_BITS * oTD.n_samples_TAG_BIT;
end

% compute maximum amplitude preamble
n_samples_preamble = (oTD.TAG_PREAMBLE_BITS * oTD.n_samples_TAG_BIT + oTD.n_samples_TAG_BIT/2);
ini_frame = ini_bits - n_samples_preamble;

% Use preamble also to compute power
magn_squared_samples = abs(in).^2;
maximum = max(magn_squared_samples(ini_frame : 1 : ...
                        ini_frame + n_samples_preamble + n_samples_bits));

% Compute average for samples that exceed threshold * max(RN16_raw)
counter = 0;
for i = (ini_frame:1:ini_frame+n_samples_bits)
    if (magn_squared_samples(i) > maximum*threshold)
        average = average + magn_squared_samples(i);
        counter = counter + 1;
    end
end
average = average / counter;

% Correct for Matched filter power increase (assume 20-25 taps)
average = average / 20;

power_dbm = 10*log10(average^2 / Zc) + 30;
