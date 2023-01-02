clc
clear
close        % figure close

p = [1 0 0 1]

code = Hamming_Enc(p')
modulated_signal = Two_PAM_mod(code, 16)
sigma_v = 2

received_signal = AWGN_channel(modulated_signal, sigma_v)

demodulated_signal = Two_PAM_dem(received_signal, sigma_v, 20)

estimation = Hamming_DEC(demodulated_signal) % estiamtion result




