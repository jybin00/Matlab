clc
clear
close

input = [1 0 0 1];

parity = Convolution_code_tail_bit(input);
modulated_signal = Sixteen_QAM(parity);
received_signal = AWGN_channel(modulated_signal);

