clear all;
close all
clc;

addpath('utils/');

SNR = 30;

Niter = 50;
N_list = (1e3:1e3:5e4);
% fo_list = [10 (50:50:200)];  % selected freq offset in Hz
fo_list = (1:2:8);  % selected freq offset in Hz

per_old = Inf;

CFO_est_list = zeros(numel(fo_list),numel(N_list));
for fo = fo_list
    for N = N_list
        CFO_est_temp = zeros(1,Niter);
        parfor iter = 1:Niter
            po = 2*pi*rand(1);  % phase offset
            CFO_est_temp(iter) = test_CFOEst_Theory(N,SNR,fo,po);
        end
        CFO_est_list(fo==fo_list,N==N_list) = mean(CFO_est_temp);
        fprintf('CFO: %d , N: %d - CFO error: %.4f Hz\n',...
                 fo,N,CFO_est_list(fo==fo_list,N==N_list));
        var(CFO_est_temp)/mean(CFO_est_temp);
    end
end

%% plotting
figure
loglog(N_list,CFO_est_list,'linewidth',1.5);
xlabel('N samples');
ylabel('CFO est error (Hz)');
grid on;
set(gca,'FontWeight','bold','fontSize',12);