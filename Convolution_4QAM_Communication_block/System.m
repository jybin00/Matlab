clc
clear
close

input = [1 0 0 1];

parity = Convolution_code(input);
modulated_signal = Four_QAM(parity);
received_signal = AWGN_channel(modulated_signal);

