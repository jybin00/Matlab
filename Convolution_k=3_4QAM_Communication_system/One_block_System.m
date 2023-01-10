clc
clear
close
tic
N_m_bit = 50;                                                      % Number of message bit
N_frame = 800;                                                 % Number of frame
test_bit = randsrc(N_frame, N_m_bit, [0 1]);       % test bit generation



Eb_db_final = 22;
Eb_No_db_sim = zeros(1, Eb_db_final);
sigma_v = 2;
FER = zeros(1, Eb_db_final);
error = zeros(1, Eb_db_final);
modulated_output = zeros(1, 2*(N_m_bit+2));


parfor Eb_db = 2: Eb_db_final             % Eb를 바꿔가면서 계산
    for j = 1 : N_frame                       % frame 개수만큼 계산

        % 1~ message bit || 1+ message bit ~ 2* message bit ....
        % encoded input = (message_bit + tail bit) * 2
        encoded_input = Convolution_code(test_bit(j, :));   

        % 4QAM => 2 signal -> 1 symbol
        for i = 1: (N_m_bit+2)
            modulated_output(1, 2*(i-1)+1 : 2*i) = four_QAM( encoded_input(1, 2*(i-1)+1 : 2*i), 10^(Eb_db/10) );
        end

        % Signal transmitt through AWGN channel with noise variance sigma_v
        received_signal = AWGN_Channel(modulated_output, sigma_v);

        % ML demodulation
        for i = 1: (N_m_bit+2)
            demodulated_output(1, 2*(i-1)+1 : 2*i)=Demodulation(received_signal(1, 2*(i-1)+1 : 2*i));
        end

        decoding = Viterbi_decoding(demodulated_output, N_m_bit);
        a = nnz(input-decoding);
        if a > 0
            FER(1, Eb_db)= FER(1, Eb_db) + 1;
        end
        error(1, Eb_db) = error(1, Eb_db) + a;
        
    end
    Eb_No_db_sim(1, Eb_db) = Eb_db - 10*log10(2*sigma_v^2);
end
BER = error / length(test_bit);
FER = FER / N_frame;
toc

%%
Eb_of_No_db = -1:0.1:15;
%figure(1), semilogy(Eb_of_No_db, 4*(1-1/2)*qfunc(sqrt(10.^(Eb_of_No_db/10) ) ), 'r--' );
%hold on
xlabel("Eb/No"); 
ylabel('BER');
%plot(Eb_No_db_sim, BER, 'bo-')
%hold on
semilogy(Eb_of_No_db, 1- (1-4*(1-1/2)*qfunc(sqrt(10.^(Eb_of_No_db/10) ) )).^(N_m_bit), 'r--' );
hold on
plot(Eb_No_db_sim, FER, 'bo-')

axis([-1 15 0.5*10^-6 1])
xticks([-1:2:15])
grid on

legend('Uncoded 4QAM BER', 'Hard Viterbi v = 2')
legend('Uncoded 4QAM FER', 'Hard Viterbi v = 2')