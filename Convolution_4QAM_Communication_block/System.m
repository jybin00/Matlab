clc
clear
close

input = randsrc(1, 96, [0 1]);
%input =   [0  1  0  0  0  1  0  0  1  1  1  0  0  1  1  1  0  0  0  1  1  1  0  0  0  0 1 1 1 1  1 1  1  0  1  0  1  0  0  1 0  1 1  0  1 0  1  1  0  0  0  0 0 1 0 1 0  1 0 1 1 1 0 1 1 1 0 0 1 0 0 0  1 1 0 1 1 1 1 1 1 0 1 0 1 1 0 1 0 0  0  1 1 0  1  1]
encoded_input = Convolution_code_tail_bit(input, 96);
for i = 1: length(encoded_input)/4
    modulated_output(1, 1+2*(i-1):2*i) = Sixteen_QAM(encoded_input(1, 1+4*(i-1):4*i), 100);
end
Channel_output = AWGN_channel(modulated_output, 2);
for i = 1:length(Channel_output)/2
    demodulated_output(1, 1+4*(i-1):4*i) = ML_demodulation(Channel_output(1,1+2*(i-1):2*i), 100);
end
decoding = Viterbi_decoding(demodulated_output, 96);

nnz(input-decoding)
