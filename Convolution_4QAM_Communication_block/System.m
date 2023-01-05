clc
clear
close

input = randsrc(1, 96, [0 1]);
%input =   [0  1  0  0  0  1  0  0  1  1  1  0  0  1  1  1  0  0  0  1  1  1  0  0  0  0 1 1 1 1  1 1  1  0  1  0  1  0  0  1 0  1 1  0  1 0  1  1  0  0  0  0 0 1 0 1 0  1 0 1 1 1 0 1 1 1 0 0 1 0 0 0  1 1 0 1 1 1 1 1 1 0 1 0 1 1 0 1 0 0  0  1 1 0  1  1]
encoded_input = Convolution_code_tail_bit(input, 96);
decoding = Viterbi_decoding(encoded_input, 96);

nnz(input-decoding)
