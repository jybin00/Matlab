clc
clear
close all

%% System parameter setting, memory allocation
M = 2;                                  % Modulation order (1: BPSK, 2: QPSK)
N_bit = 1000000;                % Number of bit
N_sym = N_bit/M;                % Number of symbol
SNR_dB = 0:2:20;                 % SNR (dB scale)
SNR = 10.^(SNR_dB/10);      % SNR (linear scale)
Noise_var = 1./SNR;             % Noise variance (linear scale)
h_Multipath = [1 0.4 0.3 0.2 0.1];  % Multi-path tap gain

err_count_Rayleigh = zeros(1,length(SNR_dB)); 
err_count_Multipath = zeros(1,length(SNR_dB));

%% Main loop
for n=1:length(SNR)           % SNR loop
    x_i = randi(2,[1,N_bit])-1;
    x_k = Modulation(M,N_sym,x_i);
    n_k = sqrt(Noise_var(n))*(1/sqrt(2))*(randn(1,N_sym)+randn(1,N_sym)*1i);
    y_k_Rayleigh = Rayleigh_Fading(N_sym,x_k)+n_k; 
    y_k_Multipath = Multipath_Fading(N_sym,x_k,h_Multipath)+n_k;
    x_i_hat_Rayleigh = Demodulation(M,N_sym,y_k_Rayleigh); x_i_hat_Multipath = Demodulation(M,N_sym,y_k_Multipath);
 for i=1:N_bit
     if x_i_hat_Rayleigh(i) ~= x_i(i)
         err_count_Rayleigh(1,n) = err_count_Rayleigh(1,n)+1;
     end
     if x_i_hat_Multipath(i) ~= x_i(i)
            err_count_Multipath(1,n) = err_count_Multipath(1,n)+1;
     end
  end
  err_count_Rayleigh(1,n) = err_count_Rayleigh(1,n)/N_bit;
  err_count_Multipath(1,n) = err_count_Multipath(1,n)/N_bit;
end

%% Result
if M == 1
    [BER_BPSK, SER_BPSK] = berfading(SNR_dB,'psk',2,1);
elseif M == 2
    [BER_QPSK, SER_QPSK] = berfading(SNR_dB-3,'psk',4,1);
end
figure(1)
if M == 1
    semilogy(SNR_dB,BER_BPSK,'bo') 
elseif M == 2
    semilogy(SNR_dB,BER_QPSK,'bo') 
end
hold on
grid on
semilogy(SNR_dB,err_count_Rayleigh,'b')
hold on
semilogy(SNR_dB,err_count_Multipath,'r')
axis([0 20 10^-5 1])
ylabel('BER','fontsize',12,'fontname','Times New Roman') 
xlabel('SNR[dB]','fontsize',12,'fontname','Times New Roman') 
legend('Rayleigh-analysis','Rayleigh-simulation','Multipath-simulation')





