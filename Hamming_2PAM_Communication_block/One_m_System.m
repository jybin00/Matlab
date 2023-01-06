clc
clear
close        % figure close

p = randsrc(1, 4, [1 0])

code = Hamming_Enc(p')
modulated_signal = Two_PAM_mod(code, 16)
sigma_v = 2

received_signal = AWGN_channel(modulated_signal, sigma_v)

%demodulated_signal = Two_PAM_dem(received_signal, sigma_v, 20)

%estimation = Hamming_DEC(demodulated_signal) % estiamtion result

estimation = Soft_decision_DEC(received_signal', 16)
p




