function output = f_CFOEst_alt(input, oMF, fs, CFO)

%% Retrieve the portion that we wish to use to estimate the CFO
my_init = 4300;  % hardcoded for sample IQ file
my_end = 6098;  % hardcoded for sample IQ file
range_CW = (my_init:my_end);
input_CFOEst = input(range_CW);
N = numel(input);
s_CFO = exp(1i*2*pi*(0:N-1)*CFO);

%% Estimate the CFO
x = (1:numel(range_CW))./fs;
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
CFO_est_Re = fo_est_Re/fs;
CFO_est_Im = fo_est_Im/fs;
phase_delta = 2*pi*CFO_est_Re*(my_init-1);
s_CFO_corr = offset_Re + (ampl_Re/oMF.n_taps) .* sin(2*pi*(0:N-1)*CFO_est_Re + phase_delta) + ...
             offset_Im + (ampl_Im/oMF.n_taps) .* 1i.*sin(2*pi*(0:N-1)*CFO_est_Im + phase_delta);

%% Compare CFO sequences
s_CFO_decim = decimate([s_CFO zeros(1,oMF.n_taps)],oMF.decim);
figure;
subplot(2,1,1); hold on;
plot((1:numel(s_CFO_decim))./fs.*1e3,real(s_CFO_decim),'linewidth',1.5);
plot((1:numel(s_CFO_corr))./fs.*1e3,real(s_CFO_corr),'linewidth',1.5);
ylabel('Real (amplitude)')
xlabel('Time (ms)');
legend('original CFO seq','estimated CFO seq');
set(gca,'FontWeight','bold','fontSize',9);
subplot(2,1,2); hold on;
plot((1:numel(s_CFO_decim))./fs.*1e3,imag(s_CFO_decim),'linewidth',1.5);
plot((1:numel(s_CFO_corr))./fs.*1e3,imag(s_CFO_corr),'linewidth',1.5);
ylabel('Imaginary (amplitude)')
xlabel('Time (ms)');
legend('original CFO seq','estimated CFO seq');
set(gca,'FontWeight','bold','fontSize',9);        
         
output = input.*conj(s_CFO_corr);