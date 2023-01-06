clc         % command line clear
clear      % variable clear
close     % figure close


N_frame = 5000;                                   % frame 수
input = randsrc(1, 96*N_frame, [0 1]);       % input의 개수

encoded_input = Convolution_code_tail_bit(input, 96);
for i = 1: length(encoded_input)/4
    modulated_output(1, 1+2*(i-1):2*i) = Sixteen_QAM(encoded_input(1, 1+4*(i-1):4*i), 100);
end
Channel_output = AWGN_channel(modulated_output, 2);
for i = 1:length(Channel_output)/2
    demodulated_output(1, 1+4*(i-1):4*i) = ML_demodulation(Channel_output(1,1+2*(i-1):2*i), 100);
end
for i = 1:length(demodulated_output)/212
    decoding(1, 1+96*(i-1):96*i) = Viterbi_decoding(demodulated_output(1, 1+212*(i-1):212*i), 96);
end

%nnz(input-decoding)
