clc
clear
tic
% 매트랩 내장 함수인 vitdec을 이용해서 실제로 제대로 구현한 것인지 확인해보기.

%system parameter
Eb_No_dB = (0:15)';
SNR_dB = Eb_No_dB;
BER_hard = zeros(length(Eb_No_dB));
BER_soft = zeros(length(Eb_No_dB));

% input
number_of_msg = 10e5;
msg = randi([0 1], 1, number_of_msg);

% convolution code
trellis = poly2trellis(3, [7 5]);
coded_bit = convenc(msg, trellis);

% mapping
symbol = 2*coded_bit - 1;

for i = 1: length(Eb_No_dB)
    % channel
    received_bit = awgn(symbol, SNR_dB(i));
    
    % demapping
    demodulated_output = received_bit > 0;
    
    % decoding
    decoding_hard = vitdec(demodulated_output, trellis, 50, 'trunc', 'hard');
    decoding_soft =  vitdec(received_bit,               trellis, 50, 'trunc', 'unquant');
    
    % FER, BER measure
    BER_hard(i) = nnz(decoding_hard - msg) / number_of_msg;
    BER_soft(i) = nnz(decoding_soft - msg) / number_of_msg;
end
toc

%%
close all
Eb_of_No_db = -1:0.1:12;
% theorical uncoded BER
semilogy(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), 'r--' );
hold on
grid on

% sim BER
plot(Eb_No_dB, BER_hard, 'kx-')

plot(Eb_No_dB, BER_soft, 'bo-')

axis([0 12 10^-5 1])
xticks(0:2:12)

xlabel("Eb/No"); 
ylabel('BER');

%legend('Uncoded 2PAM BER', 'Hard wrong Viterbi v = 2, m = 50', 'Hard correct Viterbi v = 2, m = 50')
legend('Uncoded 2PAM BER', 'Hard Viterbi v = 2', 'Soft Viterbi v = 2')

