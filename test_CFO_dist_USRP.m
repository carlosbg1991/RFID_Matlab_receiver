clear all; close all; clc;

addpath('utils/');

fileName = '~/git/rfid-mimo/gr-rfid/misc/data/usrp_sampl_w_cfo.out';
% fileName = '~/git/rfid-mimo/gr-rfid/misc/data/usrp_sampl_w_cfo_2.out';
% fileName = '~/git/rfid-mimo/gr-rfid/misc/data/usrp_sampl_w_cfo_3.out';
FS = 2e6;
fid = fopen(fileName, 'r');
data = fread(fid,inf,'single');
raw_IQ = (data(1:2:end) + 1i.*data(2:2:end)).';

% figure; hold on;
% plot((1:numel(raw_IQ)).*(1e3/FS),real(raw_IQ));
% plot((1:numel(raw_IQ)).*(1e3/FS),imag(raw_IQ));
% xlabel('Time (ms)');
% set(gca,'FontWeight','bold','fontSize',10);

CFO_exp = 150;
CFO_max = 500;  % in Hz
Y = real(raw_IQ);  % input only expected N number of cycles
X = (1:numel(Y))./FS;

N_sampl_bulk = 4*1.5*(FS/CFO_exp);  % the FFT needs 1.5 times the samples per
                                    % cycle. We double it in case the CFO
                                    % is around twice the expected.
N_bulks = ceil(numel(Y)/N_sampl_bulk);

f_peaks_list = zeros(1,N_bulks);
for b = 1:N_bulks
    idx_range = (1+N_sampl_bulk*(b-1):N_sampl_bulk*b);
    x = X(idx_range);
    y = Y(idx_range);
    f_peaks_list(b) = coarse_freq_est(x,y);
    if mod(b,10)==0
        fprintf('%.3f\n',100*b/N_bulks);
    end
end
figure; hold on;
histogram(f_peaks_list,'BinWidth',0.5);

function fPeak = coarse_freq_est(x,y)
    pi2=2*pi;
    NumSamples=length(x);
    T=x(2)-x(1);
    offs=mean(y);%DC value
    y_m=y-offs;%FFT much better without offset
    n = 128*2^nextpow2(NumSamples);%heavy zero padding
    Y = fft(y_m,n);%Y(f)
    n2=floor(n/2);
    P2 = abs(Y/NumSamples);
    P1 = P2(1:n2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    fs = (0:n2)/n/T;% frequency scale
    % %FFT parameters at peak
    [maxFFT,maxFFTindx]=max(P1);%Peak magnitude and location
    fPeak=fs(maxFFTindx);% f at peak
end