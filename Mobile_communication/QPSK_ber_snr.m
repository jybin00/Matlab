clc
clear all
close all

%% System parameter setting, memory allocation
N_bit = 3000000;
N_sym = N_bit/2;
EbNo_db = -1:1:9;
SNR_db = EbNo_db + log2(2);
SNR = 10.^(SNR_db/10);
N0 = 2./SNR;
err_count_AWGN = zeros(1,length(SNR_db));

%% Main loop
for n=1:length(SNR) % SNR loop
    x_i = randi(2,[1,N_bit])-1; %% Binary bit generation 
    x_k = QPSK_Modulation2(N_sym,x_i);

    noise = sqrt(N0(n)/2)*(1/sqrt(2))* (randn(1,N_sym)+ randn(1,N_sym)*1i);  %% Additive White Gaussian Noise

    y_AWGN = x_k + noise; %% Received signal 
    x_i_hat_AWGN= QPSK_Demodulation2(N_sym,y_AWGN); 
    %% BER Count for QPSK over AWGN
    err_count_AWGN(1,n) = sum(mod(x_i_hat_AWGN ~= x_i,2));
    err_count_AWGN(1,n) = err_count_AWGN(1,n)/N_bit;
end

%% Result
figure(1)
theoryQPSK = 1/2*(erfc(sqrt(SNR)/sqrt(2))); % theoretical BER 
semilogy(SNR_db, theoryQPSK, 'k')
hold on
grid on
semilogy(SNR_db, err_count_AWGN, 'bo')
hold on
axis([0 10 10^-4 2*10^-1]) 
ylabel('BER','fontsize',12,'fontname','Times New Roman') 
xlabel('SNR[dB]','fontsize',12,'fontname','Times New Roman') 
legend('Theoretical','AWGN Simulation')
