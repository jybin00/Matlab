clc         % command line clear
clear      % variable clear
close     % figure close
tic

N_frames = 4000000;                             % number of symbols
E_bit_db = 21;                                          % bit energy [dB]
b_error = zeros(E_bit_db, 1);                    % error bit array
f_error = zeros(E_bit_db, 1);                     % error bit array
N_f_sim(1 : E_bit_db, :) = N_frames;
Eb_of_No_dB = zeros(E_bit_db, 1);            % Eb of No[dB] aaray

sigma_v = 2;                                             % noise variance  

generated_sequnece = randsrc(N_frames, 4, [0 1]);         % #of frame X message bits

parfor Eb_db = [8 : E_bit_db]                            % Eb of No range
    for i = 1 : N_frames

        input = generated_sequnece(i, :);                                                                                   % seperate 4bits
        code = Hamming_Enc(input');                                                                                         % hamming encoding 
        modulated_signal = Two_PAM_mod(code, (4/7)*10^(Eb_db/10));                                   % 2PAM modulation
        received_signal = AWGN_channel(modulated_signal, sigma_v);                                       % AWGN channel 
        % demodulated_signal = Two_PAM_dem(received_signal, sigma_v, 10^(Eb_db/10));       % LLR demodulation
        % estimation = Hamming_DEC(demodulated_signal);                                                      % DEC
        %-----------------------------------------------------------------------------------
        estimation = Soft_decision_DEC(received_signal', (4/7)*10^(Eb_db/10));                            % soft decision

        bit_error = nnz(input - estimation);
        b_error(Eb_db,1)  = b_error(Eb_db,1)   + bit_error;      % count error bit
        if bit_error > 0
            f_error(Eb_db,1) = f_error(Eb_db,1)  + 1;                 % count error frame 
        end
        if b_error(Eb_db,1) > 200
            N_f_sim(Eb_db, 1) = i;
            disp(N_f_sim(Eb_db, 1));
            break;
        end
    end
    BER(Eb_db, 1) = b_error(Eb_db,1) / (N_f_sim(Eb_db, 1) * 4);       % BER calculation alog Eb
    FER(Eb_db, 1) = f_error(Eb_db,1)  / (N_f_sim(Eb_db, 1));            % FER calculation along Eb
    Eb_of_No_dB(Eb_db) = Eb_db - 10*log10(2*(sigma_v^2));         % Eb[db] - No[db] = Eb/No [db]
end

toc

%% drawing graph
close
Eb_No = -2:0.01:12;
figure(1),semilogy( Eb_No, (qfunc ( sqrt(2*10.^((Eb_No)/10)) )), 'b--' );   % uncoded 2PAM BER graph
hold on

semilogy( Eb_No, 1- (1-qfunc ( sqrt(2*10.^(Eb_No/10)))).^4, 'g--' );      % uncoded 2PAM FER graph
xlabel('Eb/No [dB]'), ylabel('FER');
hold on

plot(Eb_of_No_dB, BER, '-kx');                      % BER
hold on

plot(Eb_of_No_dB, FER, '-ro')                       % FER
hold on
grid,axis([-0.2 12 0.5*10^(-6) 1]); 
legend("Uncoded 2PAM BER", "Uncoded 2PAM FER", "Simulation BER", "Simulation FER");





