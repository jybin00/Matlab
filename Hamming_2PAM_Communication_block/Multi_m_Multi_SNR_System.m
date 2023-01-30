clc         % command line clear
clear      % variable clear
close     % figure close
tic

N_frames = 5000000;                             % number of symbols
EbNo_db = (0:10)';                                        % bit energy [dB]
b_error = zeros(length(EbNo_db), 1);                   % error bit array
b_error_h = zeros(length(EbNo_db), 1);                   % error bit array

f_error = zeros(length(EbNo_db), 1);                    % error bit array
f_error_h = zeros(length(EbNo_db), 1);                    % error bit array
N_f_sim(1 : length(EbNo_db), :) = N_frames;

k = 1;
rate = 4/7;


generated_sequnece = randsrc(N_frames, 4, [0 1]);         % #of frame X message bits

parfor n = 1:length(EbNo_db)                                         % Eb of No range
    snrdB = EbNo_db(n) + 10*log10(k*rate);
    noise_Var = 10.^(-snrdB/10);
    sigma = sqrt(noise_Var/2)
    for i = 1 : N_frames

        input = generated_sequnece(i, :);                                                                        % seperate 4bits
        code = Hamming_Enc(input');                                                                              % hamming encoding 
        modulated_signal = Two_PAM_mod(code, 1);                                                      % 2PAM modulation

        received_signal = AWGN_channel(modulated_signal, sigma);                            % AWGN channel
        %received_signal = awgn(modulated_signal, snrdB, 'measured');

        demodulated_signal = Two_PAM_dem(received_signal);                                      % LLR demodulation
        estimation_h = Hamming_DEC(demodulated_signal);                                          % DEC
        %-----------------------------------------------------------------------------------
        estimation_s = Hamming_DEC(Soft_decision_DEC(received_signal', 1));               % soft decision

        bit_error_s = nnz(input - estimation_s);
        bit_error_h = nnz(input - estimation_h);

        b_error(n,1)  = b_error(n,1)   + bit_error_s;             % count error bit
        b_error_h(n,1)  = b_error_h(n,1)   + bit_error_h;      % count error bit
        if bit_error_s > 0
            f_error(n,1) = f_error(n,1)  + 1;                         % count error frame
        end
        if bit_error_h > 0
            f_error_h(n,1) = f_error_h(n,1)  + 1;
        end
        if b_error(n,1) > 1000
            N_f_sim(n, 1) = i;
            disp(N_f_sim(n, 1));
            break;
        end
    end
    BER(n, 1) = b_error(n,1) / (N_f_sim(n, 1) * 4);            % BER calculation alog Eb
    BER_h(n, 1) = b_error_h(n,1) / (N_f_sim(n, 1) * 4);     % BER calculation alog Eb
    FER(n, 1) = f_error(n,1)  / (N_f_sim(n, 1));                 % FER calculation along Eb
    FER_h(n, 1) = f_error_h(n,1)  / (N_f_sim(n, 1));   
end
toc

%% drawing graph
close
Eb_No_db = -2:0.1:12;
figure(1),semilogy( Eb_No_db, (qfunc ( sqrt(2*10.^((Eb_No_db)/10)) )), 'b--' );   % uncoded 2PAM BER graph
hold on
k = 1;
theorical_fer = zeros(length(Eb_No_db), 1);
for i = Eb_No_db
    p = qfunc(sqrt((8/7)*10^(i/10)));
    theorical_fer(k,1) = (1- ( (1-p)^7 + 7*p*(1-p)^6 ))*(3/7);
    k = k + 1;
end
semilogy(Eb_No_db, theorical_fer)
hold on 

plot(EbNo_db, BER_h, 'ko-');                      % hard BER
hold on

%plot(EbNo_db, BER, '-kx');                          % soft BER

grid,axis([-0.2 12 0.5*10^(-6) 1]); 
legend("Uncoded 2PAM BER", "Theoretical Hamming BER", "Hard decision BER", "Soft decision BER");

%%
Eb_No_db = -2:0.01:12;
semilogy( Eb_No_db, 1- (1-qfunc ( sqrt(2*10.^(Eb_No_db/10)))).^4, 'g--' );      % uncoded 2PAM FER graph
xlabel('Eb/No [dB]'), ylabel('FER');
hold on

plot(EbNo_db, FER_h, '-ro')                       % FER

grid,axis([-0.2 12 0.5*10^(-6) 1]); 
legend("Uncoded 2PAM FER", "Simulation FER");
