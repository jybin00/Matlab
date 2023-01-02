clc
clear
close
syms y

N_frames = 100000;  % number of symbols
correct = zeros(2, 1);
error = zeros(2, 1);
c_or_e = zeros(2, 1);
c_or_e2 = zeros(2, 1);
generated_sequnece = randsrc(4*N_frames, 1, [0 1]);

E_b = 10;
sigma_v = 2;             % noise variance 

for i = 1 : N_frames
    input = generated_sequnece(4*(i-1)+1 : 4*i)';

    code = Hamming_Enc(input');

    modulated_signal = Two_PAM_mod(code, E_b);

    received_signal = AWGN_channel(modulated_signal, sigma_v);
    
    demodulated_signal = Two_PAM_dem(received_signal, sigma_v, E_b);
    
    estimation = Hamming_DEC(demodulated_signal);

    c_or_e = BER_analysis(input, estimation);
    correct(1) = correct(1) + c_or_e(1);
    error(1)    = error (1)   + c_or_e(2);

    c_or_e2 = FER_analysis(input, estimation);
    correct(2) = correct(2) + c_or_e2(1);
    error(2)   = error(2)    + c_or_e2(2);
end
BER = error(1)/(correct(1) + error(1))
FER = error(2)/(correct(2) + error(2))
Eb_of_No_dB = 10*log10(E_b/(2*(sigma_v^2)))


%% drawing graph
Eb_No = -2:0.01:14;
figure(1),semilogy( Eb_No, (qfunc ( sqrt(2*10.^((Eb_No)/10)) )), 'r--' );  % uncoded 2PAM graph
xlabel('Eb/No [dB]'), ylabel('BER')
hold on
scatter(solve(Shannon_Cap(y, sigma_v) == 4/7, y), [10^(-7) 1])
hold on
scatter(Eb_of_No_dB, FER, 20, 'magenta', "X");
hold on
scatter(Eb_of_No_dB, BER, 20, 'green', 'filled');
grid,axis([-2 14 10^(-7) 1]);
legend('Uncoded 2PAM result', 'Shannon limits','Simulation BER', 'Simulation FER');





