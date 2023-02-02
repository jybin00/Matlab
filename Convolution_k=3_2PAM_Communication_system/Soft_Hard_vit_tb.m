clc
clear
tic
N_m_bit = 50;                                                    % Number of message bit
N_frame = 400000;                                            % Number of frame
test_bit = randsrc(N_frame, N_m_bit, [0 1]);       % test bit generation


Eb_db_final = 19;
Eb_No_db_sim = zeros(1, Eb_db_final);

FER_soft = zeros(1, Eb_db_final);
FER_hard = zeros(1, Eb_db_final);

BER_soft = zeros(1, Eb_db_final);
BER_hard = zeros(1, Eb_db_final);

error_soft = zeros(1, Eb_db_final);
error_hard = zeros(1, Eb_db_final);

N_f_sim(1,1:Eb_db_final) = N_frame;
sigma_v = 2;
demodulated_output = zeros(1, 2*(N_m_bit+2));


parfor Eb_db = 8: Eb_db_final          % Eb를 바꿔가면서 계산
    for j = 1 : N_frame                       % frame 개수만큼 계산

        % 1~ message bit || 1+ message bit ~ 2* message bit ....
        % encoded input = (message_bit + tail bit) * 2
        input = test_bit(j, :);
        encoded_input = Convolution_code(input);  
        demodulated_output = zeros(1, 2*(N_m_bit+2));

        encoded_input = reshape(encoded_input, [2 N_m_bit+2]);
        encoded_input = encoded_input';
        modulated_output = zeros(1, 2*(N_m_bit+2));

        % modulation -> channel -> demodulation
        for i = 1: (N_m_bit+2)
            % 4QAM => 2 signal -> 1 symbol
            modulated_output(1, i) = four_QAM( encoded_input(i, :), 10^(Eb_db/10) );
        end

        % Signal transmitt through AWGN channel with noise variance sigma_v
        received_signal = AWGN_Channel(modulated_output, (sigma_v/sqrt(2)));

        % ML demodulation
        for i = 1: (N_m_bit+2)
            demodulated_output(1, 2*(i-1)+1 : 2*i)=Demodulation(received_signal(1,i));
        end

        received_signal= [real(received_signal); imag(received_signal)];
        decoding_soft = Wrong_Viterbi_soft_decoding(received_signal, N_m_bit, 10^(Eb_db/10));
        decoding_hard = Wrong_Viterbi_decoding(demodulated_output, N_m_bit);

        error_s = nnz(input-decoding_soft);
        error_h = nnz(input-decoding_hard);
        if error_s > 0
            FER_soft(1, Eb_db)= FER_soft(1, Eb_db) + 1;
            error_soft(1, Eb_db) = error_soft(1, Eb_db) + error_s;
        end
        if error_h > 0
            FER_hard(1, Eb_db)= FER_hard(1, Eb_db) + 1;
            error_hard(1, Eb_db) = error_hard(1, Eb_db) + error_h;
        end

        if error_soft(1, Eb_db) >400
            N_f_sim(1, Eb_db) = j;
            break
        end
        
    end
    Eb_No_db_sim(1, Eb_db) = Eb_db - 10*log10(2*sigma_v^2);
    BER_soft(1, Eb_db) = error_soft(1, Eb_db) / (N_f_sim(1, Eb_db) * N_m_bit);
    BER_hard(1, Eb_db) = error_hard(1, Eb_db) / (N_f_sim(1, Eb_db) * N_m_bit);

    FER_soft(1, Eb_db) = FER_soft(1, Eb_db) / N_f_sim(1, Eb_db);
    FER_hard(1, Eb_db) = FER_hard(1, Eb_db) / N_f_sim(1, Eb_db);
end
% BER = error / numel(test_bit);
% FER = FER / N_frame;
toc
%%
close all
Eb_of_No_db = -1:0.1:15;
% theorical BER
figure(1), semilogy(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), 'r--' );
hold on

% sim BER
plot(Eb_No_db_sim, BER_hard, 'kx-')
hold on
plot(Eb_No_db_sim, BER_soft, 'bo-')

axis([0 14 0.5*10^-6 1])
xticks(0:2:12)
grid on

xlabel("Eb/No"); 
ylabel('BER');

legend('Uncoded 4QAM BER', 'Hard Viterbi v = 2, m = 50', 'Soft Viterbi v = 2, m = 50')
%legend('Uncoded 4QAM BER', 'Hard Viterbi v = 2, m = 25', 'Soft Viterbi v = 2, m = 25')
legend('Uncoded 4QAM BER', 'Soft Viterbi v = 2, m = 50')

%%
close all
Eb_of_No_db = -1:0.1:15;
figure(2)
% theorical FER
semilogy(Eb_of_No_db, 1- (1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) )).^(50), 'r--' );
hold on

% sim FER
%plot(Eb_No_db_sim, FER_hard, 'bo-')
%hold on
plot(Eb_No_db_sim, FER_soft, 'bo-')

axis([0 14 10^-5 10])
xticks([0:2:12])
grid on

xlabel("Eb/No"); 
ylabel('FER');

%legend('Uncoded 4QAM FER', 'Hard Viterbi v = 2, m = 50', 'Soft Viterbi v = 2, m = 50')
legend('Uncoded 4QAM FER', 'Soft Viterbi v = 2, m = 50')