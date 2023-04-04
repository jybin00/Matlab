clc
clear all
close all

%% System parameter setting, memory allocation
N_bit = 1000000;
N_sym = N_bit/2;
EbN0_dB = 0:1:8;
EbN0 = 10.^(EbN0_dB/10);
N0 = 1./EbN0;
err_count_AWGN = zeros(1,length(EbN0_dB));

%% Main loop
for n=1:length(EbN0) % SNR loop
    x_i = randi(2,[1,N_bit])-1; %% Binary bit generation 
    x_k = QPSK_Modulation(N_sym,x_i);

    noise = sqrt(N0(n)/2)*(1/sqrt(2))* (randn(1,N_sym)+ randn(1,N_sym)*1i);  %% Additive White Gaussian Noise

    y_AWGN = x_k + noise; %% Received signal 
    x_i_hat_AWGN= QPSK_Demodulation(N_sym,y_AWGN); 
    %% BER Count for QPSK over AWGN
    for i=1:N_bit
        if x_i_hat_AWGN(i) ~= x_i(i)
            err_count_AWGN(1,n) = err_count_AWGN(1,n)+1;
        end 
    end
    err_count_AWGN(1,n) = err_count_AWGN(1,n)/N_bit;
end

%% Result
figure(1)
theoryQPSK = 0.5*(erfc(sqrt(EbN0))); % theoretical BER semilogy(EbN0_dB,theoryQPSK,'k')
hold on
grid on
semilogy(EbN0_dB,err_count_AWGN,'bo')
hold on
axis([0 8 10^-4 10^-1]) 
ylabel('BER','fontsize',12,'fontname','Times New Roman') 
xlabel('EbN0[dB]','fontsize',12,'fontname','Times New Roman') 
legend('Theoretical','AWGN Simulation')