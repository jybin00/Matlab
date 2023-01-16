% soft decision
clc
clear
close
%input = randsrc(1, 10, [0 1])
input = [1 0 1 1 0 0 1 1 1 1 0 0]
encoded_input = Convolution_code(input)
for i = 1:length(encoded_input)/2
    modulated_output(1,i)=four_QAM(encoded_input(1,2*(i-1)+1 : 2*i), 20);
end
received_signal = AWGN_Channel(modulated_output,2)
received_signal=[real(received_signal);imag(received_signal)]
decoding = Viterbi_soft_decoding(received_signal, 12, 20)
a = nnz(input-decoding)