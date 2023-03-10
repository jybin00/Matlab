clc
clear
tic
N_m_bit = 25;                                                    % Number of message bit
N_frame = 40000;                                            % Number of frame
test_bit = randsrc(N_frame, N_m_bit, [0 1]);       % test bit generation


Eb_db_final = 19;
Eb_No_db_sim = zeros(1, Eb_db_final);
sigma_v = 2;
FER = zeros(1, Eb_db_final);
BER = zeros(1, Eb_db_final);
error = zeros(1, Eb_db_final);
N_f_sim(1,1:Eb_db_final) = N_frame;
demodulated_output = zeros(1, 2*(N_m_bit+2));


parfor Eb_db = 6: Eb_db_final             % Eb를 바꿔가면서 계산
    for j = 1 : N_frame                       % frame 개수만큼 계산

        % 1~ message bit || 1+ message bit ~ 2* message bit ....
        % encoded input = (message_bit + tail bit) * 2
        input = test_bit(j, :);
        encoded_input = Convolution_code(input);  
        demodulated_output = zeros(1, 2*(N_m_bit+2));

        encoded_input = reshape(encoded_input, [2 N_m_bit+2]);
        encoded_input = encoded_input';

        % modulation -> channel -> demodulation
        for i = 1: (N_m_bit+2)
            % 4QAM => 2 signal -> 1 symbol
            modulated_output = four_QAM( encoded_input(i, :), 10^(Eb_db/10) );

            % Signal transmitt through AWGN channel with noise variance sigma_v
            received_signal = AWGN_Channel(modulated_output, sigma_v/sqrt(2));

            % ML demodulation
            demodulated_output(1, 2*(i-1)+1 : 2*i)=Demodulation(received_signal);
        end

        decoding = Viterbi_decoding(demodulated_output, N_m_bit);
        a = nnz(input-decoding);
        if a > 0
            FER(1, Eb_db)= FER(1, Eb_db) + 1;
            error(1, Eb_db) = error(1, Eb_db) + a;
        end
        if error(1, Eb_db) > 200
            N_f_sim(1, Eb_db) = j;
            break
        end
    end
    Eb_No_db_sim(1, Eb_db) = Eb_db - 10*log10(2*sigma_v^2);
    BER(1, Eb_db) = error(1, Eb_db) / (N_f_sim(1, Eb_db) * N_m_bit);
    FER(1, Eb_db) = FER(1, Eb_db) / N_f_sim(1, Eb_db);
end
% FER = FER / N_frame;
toc

%%
close all
Eb_of_No_db = -1:0.1:12;
% theorical BER
figure(1), semilogy(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), 'r--' );
hold on

% label

% sim BER
plot(Eb_No_db_sim, BER, 'bo-')
%hold on
%plot(x, y, 'ko-')

% theorical FER
% semilogy(Eb_of_No_db, 1- (1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) )).^(N_m_bit), 'r--' );
% hold on
% sim FER
% plot(Eb_No_db_sim, FER, 'bo-')
 %plot(x, y, 'ko-')

axis([0 12 0.8*10^-6 1])
xticks(0:2:14)
grid on

xlabel("Eb/No"); 
ylabel('BER');

%legend('Uncoded 4QAM BER', 'Hard Viterbi v = 2, m = 50', 'Soft Viterbi v = 2, m = 50')
legend('Uncoded 4QAM BER', 'Hard Viterbi v = 2, m = 25', 'Hard Viterbi v = 2, m = 50')
%%
close all
Eb_of_No_db = 0:0.1:12;
figure(2)
% theorical FER
semilogy(Eb_of_No_db, 1- (1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) )).^(N_m_bit), 'r--' );
hold on
% sim FER
plot(Eb_No_db_sim, FER, 'bo-')
 %plot(x, y, 'ko-')

axis([0 14 10^-5 10])
xticks(0:2:14)
grid on

xlabel("Eb/No"); 
ylabel('FER');

legend('Uncoded 4QAM FER', 'Soft Viterbi v = 2, m = 50', 'Hard Viterbi v = 2, m = 50')