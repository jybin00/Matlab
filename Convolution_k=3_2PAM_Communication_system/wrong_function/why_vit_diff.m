% 왜 처음만든 비터비와 수정한 비터비의 성능이 다를까?
rng('default')
clc
clear
input = [0 1 1 1 1 0 0 1 0 1 0 0 0 1 0 0 1 0 1];

encoded_output = Convolution_code(input);
encoded_output = 2*encoded_output -1;

channel_output = awgn(encoded_output, 4, "measured");
channel_output= channel_output>0;

decoding_wrong = Wrong_Viterbi_decoding(channel_output, length(input));
decoding_correct = F_Viterbi_decoding(channel_output, length(input));

nnz(input - decoding_wrong)
nnz(input - decoding_correct)
nnz(decoding_correct - decoding_wrong)