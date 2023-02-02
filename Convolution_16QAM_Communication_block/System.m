clc         % command line clear
clear      % variable clear
close     % figure close


N_frame = 30;                                   % frame 수
N_message = 20;
input = randsrc(N_frame, N_message, [0 1]);       % input의 개수

for j = 1:N_frame
    encoded_input = Convolution_code_tail_bit(input(j , :) );
    for i = 1: length(encoded_input)/4
        modulated_output(1, i) = Sixteen_QAM(encoded_input(1, 1+4*(i-1):4*i), 50);
    end
    Channel_output = AWGN_channel(modulated_output, 1);
    
    for i = 1:length(Channel_output)
        demodulated_output(1, 1+4*(i-1):4*i) = ML_demodulation(Channel_output(1, i), 50);
    end
    
    nnz(encoded_input - demodulated_output);
    decoding = Viterbi_decoding(demodulated_output, N_message);
    
    nnz(input(j,:)-decoding)
end