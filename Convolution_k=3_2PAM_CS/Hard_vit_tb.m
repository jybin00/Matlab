clc
clear
tic
N_m_bit = 25;                                                    % Number of message bit
N_frame = 400000;                                            % Number of frame
test_bit = randsrc(N_frame, N_m_bit, [0 1]);       % test bit generation

Eb_No_db = [0:1:10]';
FER = zeros(1, length(Eb_No_db));
BER = zeros(1, length(Eb_No_db));
error = zeros(1, length(Eb_No_db));
N_f_sim(1,1:length(Eb_No_db)) = N_frame;
rate =1/2;
k = 1;

parfor i = 1: length(Eb_No_db)                    % Eb를 바꿔가면서 계산
    
    SNR_db = Eb_No_db(i) + 10*log10(k * rate);
    No = 10^(-SNR_db/10);
    sigma = sqrt(No/2);

    for j = 1 : N_frame                              % frame 개수만큼 계산
        % encoded input = (message_bit + tail bit) * 2
        input = test_bit(j, :);
        encoded_input = Convolution_code(input);  
        encoded_input = encoded_input';
        modulated_output = 2*encoded_input - 1;      % 2PAM or BPSK

        % Signal transmitt through AWGN channel with noise variance sigma_v
        received_signal = AWGN_Channel(modulated_output, sigma/sqrt(2));
        demodulated_output = zeros(1, 2*(N_m_bit+2)); 

        % demodulation
        for index = 1: 2*(N_m_bit+2)
            % ML demodulation
            demodulated_output(1, index)=Demodulation(received_signal(index));
        end

        decoding = Viterbi_decoding(demodulated_output, N_m_bit);
        a = nnz(input-decoding);
        if a > 0
            FER(1, i)= FER(1, i) + 1;
            error(1, i) = error(1, i) + a;
        end
        if error(1, i) > 1000
            N_f_sim(1, i) = j
            break
        end
    end
    BER(1, i) = error(1, i) / (N_f_sim(1, i) * N_m_bit);
    FER(1, i) = FER(1, i) / N_f_sim(1, i);
end
% FER = FER / N_frame;
toc

%%
close all
Eb_of_No_db = -1:0.1:12;
% theorical BER
figure(1), semilogy(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), 'r--' );
hold on

% sim BER
plot(Eb_No_db, BER, 'bo-')

axis([0 12 0.8*10^-6 1])
xticks(0:2:14)
grid on

% label
xlabel("Eb/No [dB]"); 
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