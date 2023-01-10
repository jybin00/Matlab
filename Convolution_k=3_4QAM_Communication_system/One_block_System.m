clc
clear
close
N_m_bit = 200;                                                      % Number of message bit
N_f_bit = 1000;                                                     % Number of frame
test_num = randsrc(1, N_m_bit * N_f_bit, [0 1]);    % test bit generation
encoded_input = Convolution_code(input);            % 


Eb_db = 22;
Eb_No_db_sim = zeros(1, Eb_db);
sigma_v = 2;


for j = 2:Eb_db
    for i = 1:length(encoded_input)/2
        modulated_output(1,2*(i-1)+1 : 2*i)=four_QAM(encoded_input(1,2*(i-1)+1 : 2*i), 10^(Eb_db/10));
    end
    received_signal = AWGN_Channel(modulated_output, sigma_v);
    for i = 1:length(modulated_output)/2
        demodulated_output(1,2*(i-1)+1 : 2*i)=Demodulation(received_signal(1,2*(i-1)+1 : 2*i));
    end
    decoding = Viterbi_decoding(demodulated_output, N_m_bit)
    a = nnz(input-decoding)
    BER = a/numel(input)
    Eb_No_db_sim(Eb_db) = Eb_db - 10*log10(2*sigma_v^2);
end



%%
Eb_of_No_db = -1:0.1:15;
figure(1), semilogy(Eb_of_No_db, 4*(1-1/2)*qfunc(sqrt(10.^(Eb_of_No_db/10) ) ), 'r--' );
xlabel("Eb/No"); 
ylabel('BER');

axis([-1 15 10^-7 1])
xticks([-1:15])

legend(['Uncoded 4QAM BER'])