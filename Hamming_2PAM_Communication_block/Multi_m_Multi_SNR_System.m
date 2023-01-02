clc         % command line clear
clear      % variable clear
close     % figure close
syms y  % calculation variable

N_frames = 50000;                                 % number of symbols
E_bit_db = 25;                                          % bit energy [dB]
correct = zeros(E_bit_db, 2);                     % correct bit array
error = zeros(E_bit_db, 2);                        % error bit array
BER = zeros(E_bit_db, 1);                          % bit error rate array
FER = zeros(E_bit_db, 1);                          % frame error rate array
Eb_of_No_dB = zeros(E_bit_db, 1);            % Eb of No[dB] aaray

sigma_v = 2;                                             % noise variance  

for Eb_db = 7 : E_bit_db                            % Eb of No range
    generated_sequnece = randsrc(1, 4*N_frames, [0 1]);         % #of frame X message bits
    for i = 1 : N_frames

        input = generated_sequnece(1, 4*(i-1)+1 : 4*i);                                                         % seperate 4bits
        code = Hamming_Enc(input');                                                                                     % hamming encoding 
        modulated_signal = Two_PAM_mod(code, 10^(Eb_db/10));                                        % 2PAM modulation
        received_signal = AWGN_channel(modulated_signal, sigma_v);                                   % AWGN channel 
        demodulated_signal = Two_PAM_dem(received_signal, sigma_v, 10^(Eb_db/10));      % LLR demodulation
        estimation = Hamming_DEC(demodulated_signal);                                                     % DEC
    
        c_or_e = BER_analysis(input, estimation);                           % correct or error bit save
        correct(Eb_db,1) = correct(Eb_db,1) + c_or_e(1);                % count correct bit
        error(Eb_db,1)    = error(Eb_db,1)    + c_or_e(2);                % count error bit

        c_or_e2 = FER_analysis(input, estimation);                         % correct or error frame save
        correct(Eb_db,2) = correct(Eb_db,2) + c_or_e2(1);              % count correct frame
        error(Eb_db,2)    = error(Eb_db,2)    + c_or_e2(2);              % count error frame
    end

    BER(Eb_db) = error(Eb_db, 1)/(correct(Eb_db, 1) + error(Eb_db, 1));     % BER calculation alog Eb
    FER(Eb_db) = error(Eb_db, 2)/(correct(Eb_db, 2) + error(Eb_db, 2));     % FER calculation along Eb
    Eb_of_No_dB(Eb_db) = Eb_db - 10*log10(2*(sigma_v^2));                   % Eb[db] - No[db] = Eb/No [db]
end


%% drawing graph
Eb_No = -2:0.01:14;
figure(1),semilogy( Eb_No, (qfunc ( sqrt(2*10.^((Eb_No)/10)) )), 'r--' );   % uncoded 2PAM BER graph
hold on

semilogy( Eb_No, 1- (1-qfunc ( sqrt(2*10.^(Eb_No/10)))).^4, 'b--' );      % uncoded 2PAM FER graph
xlabel('Eb/No [dB]'), ylabel('BER, FER');
hold on

xline(-1.6, '-.m')
hold on
    
scatter(Eb_of_No_dB, BER, 20, "green", "filled", "o");
hold on

scatter(Eb_of_No_dB, FER, 20, "red", "filled", "o");
grid,axis([-2.2 14 10^(-7) 1]), legend('Uncoded 2PAM BER', 'Uncoded 2PAM FER', 'Shannon limit', ...
    'Simulation BER', 'Simulated FER');




