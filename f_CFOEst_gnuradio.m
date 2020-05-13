function [CFO_est] = f_CFOEst_gnuradio(x,Fs,FFT_len)
% [CFO_est] = f_CFOEst_gnuradio(x,Fs,FFT_len)
% Performs tasks of coarse-grained frame synchronization by looking at the
% pulses percieved from the Query and ACK packets sent from Reader.
%
% Syntax:  [CFO_est] = f_CFOEst_gnuradio(x,Fs,FFT_len)
%
% Inputs:
%    x  - 
%    Fs - 
%    FFT_len - 
%
% Outputs:
%    CFO_est - 
%
% Example
%    N = 1e6;
%    adc_rate = 2e6;
%    fo = 150;
%    CFO = fo/adc_rate;
%    po = 0;
%    s_CFO = exp(1i*(2*pi*(0:N-1)*CFO + po));
% 
%    FFT_len = 2^15;
%    [CFO_est] = f_CFOEst_gnuradio(s_CFO,adc_rate,FFT_len);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% 
% See also: main_RX_Matlab.m

%------------- BEGIN CODE --------------
    
%% FFT
NumSamples = length(x);

offs = mean(x);%DC value
y_m = x-offs;%FFT much better without offset

n = 2^nextpow2(FFT_len);  % no zero padding
Af = Fs/FFT_len;  % frequency resolution of FFT

Y = fft(y_m, n);

n2 = floor(n/2);
P2 = abs(Y/NumSamples);
P1 = P2(1:n2+1);
P1(2:end-1) = 2*P1(2:end-1);

[~,maxFFTindx]=max(P1);%Peak magnitude and location
fs = (0:1:n2).*Af;  % frequency scale
CFO_est = fs(maxFFTindx);% f at peak

figure;  % for FFT
plot(fs,P1,'r-');
xlabel('Frequency [Hz]');
ylabel('Amplitude');


end