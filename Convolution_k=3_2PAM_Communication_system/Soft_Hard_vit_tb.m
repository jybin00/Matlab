clc
clear
tic
N_m_bit = 50;                                                 % Number of message bit
N_frame = 20000;                                             % Number of frame
test_bit = randsrc(N_frame, N_m_bit, [0 1]);     % test bit generation


Eb_No_dB = (0:10)';
Es_No_dB = Eb_No_dB + 10*log10(1/2);
FER_soft = zeros(1, length(Eb_No_dB));
FER_hard = zeros(1, length(Eb_No_dB));

BER_soft = zeros(1, length(Eb_No_dB));
BER_hard = zeros(1, length(Eb_No_dB));

error_soft = zeros(1, length(Eb_No_dB));
error_hard = zeros(1, length(Eb_No_dB));

N_f_sim(1,1:length(Eb_No_dB)) = N_frame;
demodulated_output = zeros(1, 2*(N_m_bit+2));


parfor n = 1 : length(Eb_No_dB)         % Eb를 바꿔가면서 계산
    disp(n)                                      % 진행상황 확인을 위한 인덱스
    No = 10^(-Es_No_dB(n)/10);
    sigma = sqrt(No/2);
    for j = 1 : N_frame                       % frame 개수만큼 계산

        % 1~ message bit || 1+ message bit ~ 2* message bit ....
        % encoded input = (message_bit + tail bit) * 2
        input = test_bit(j, :);
        encoded_input = Convolution_code(input);  
        demodulated_output = zeros(1, 2*(N_m_bit+2));

        encoded_input = encoded_input';
        modulated_output = zeros(2*(N_m_bit+2), 1);

        % modulation -> channel -> demodulation
        for i = 1: 2*(N_m_bit+2)
            modulated_output(i, 1) = Two_PAM(encoded_input(i), 1);
        end

        % Signal transmitt through AWGN channel with noise variance sigma_v
        received_signal = AWGN_Channel(modulated_output, sigma);

        % demodulation
        for i = 1: 2*(N_m_bit+2)
            demodulated_output(1, i)=Demodulation(received_signal(i));
        end
        
        decoding_input = reshape(received_signal, [2, length(received_signal)/2]);
        decoding_input = decoding_input';
        decoding_soft = Viterbi_soft_decoding(decoding_input, N_m_bit, 1);
        decoding_hard = Viterbi_decoding(demodulated_output, N_m_bit);

        error_s = nnz(input-decoding_soft);
        error_h = nnz(input-decoding_hard);
        if error_s > 0
            FER_soft(1, n)= FER_soft(1, n) + 1;
            error_soft(1, n) = error_soft(1, n) + error_s;
        end
        if error_h > 0
            FER_hard(1, n)= FER_hard(1, n) + 1;
            error_hard(1, n) = error_hard(1, n) + error_h;
        end

        if error_soft(1, n) > 10000
            N_f_sim(1, n) = j;
            disp(j)
            break
        end
        
    end
    BER_soft(1, n) = error_soft(1, n) / (N_f_sim(1, n) * N_m_bit);
    BER_hard(1, n) = error_hard(1, n) / (N_f_sim(1, n) * N_m_bit);

    FER_soft(1, n) = FER_soft(1, n) / N_f_sim(1, n);
    FER_hard(1, n) = FER_hard(1, n) / N_f_sim(1, n);
end
toc

%%
close all
Eb_of_No_db = -1:0.1:12;
% theorical uncoded BER
semilogy(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), 'r--' );
hold on

% sim BER
plot(Eb_No_dB, BER_hard, 'kx-')
hold on
plot(Eb_No_dB, BER_soft, 'bo-')

axis([0 12 10^-5 1])
xticks(0:2:12)
grid on

xlabel("Eb/No"); 
ylabel('BER');

%legend('Uncoded 2PAM BER', 'Hard wrong Viterbi v = 2, m = 50', 'Hard correct Viterbi v = 2, m = 50')
legend('Uncoded 2PAM BER', 'Hard Viterbi v = 2, m = 25', 'Soft Viterbi v = 2, m = 25')

%%

close all
Eb_of_No_db = -1:0.1:12;
figure(2)
% theorical FER
semilogy(Eb_of_No_db, 1- (1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) )).^(50), 'r--' );
hold on

% sim FER
plot(Eb_No_dB, FER_hard, 'kx-')
hold on
plot(Eb_No_dB, FER_soft, 'bo-')

axis([0 12 10^-5 10])
xticks([0:2:12])
grid on

xlabel("Eb/No"); 
ylabel('FER');

legend('Uncoded 2PAM FER', 'Hard Viterbi v = 2, m = 50', 'Soft Viterbi v = 2, m = 50')
% legend('Uncoded 4QAM FER', 'Soft Viterbi v = 2, m = 50')