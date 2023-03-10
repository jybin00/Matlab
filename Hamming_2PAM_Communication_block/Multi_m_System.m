clc
clear
close
syms y

N_frames = 4000000;  % number of symbols
correct = zeros(2, 1);
error = zeros(2, 1);
c_or_e = zeros(2, 1);
c_or_e2 = zeros(2, 1);
generated_sequnece = randsrc(4*N_frames, 1, [0 1]);

E_s = 20;
sigma_v = 2;             % noise variance 

for i = 1 : N_frames
    
    input = generated_sequnece(4*(i-1)+1 : 4*i)';

    code = Hamming_Enc(input');

    modulated_signal = Two_PAM_mod(code, E_s);

    received_signal = AWGN_channel(modulated_signal, sigma_v);
    
%     demodulated_signal = Two_PAM_dem(received_signal, sigma_v, E_b);
%     
%     estimation = Hamming_DEC(demodulated_signal);

    estimation = Soft_decision_DEC(received_signal', E_s);

    bit_error = nnz(input - estimation);

    error(1,1) = error(1,1) + bit_error;

    if(bit_error >0)
        error(2,1) = error(2,1) + 1;
    end

end
BER = error(1,1)/(4*N_frames);
FER = error(2,1)/N_frames;
Eb_of_No_dB = 10*log10(E_s/((4/7)*2*(sigma_v^2)));


%% drawing graph
Eb_No = -2:0.01:14;
figure(1),semilogy( Eb_No, (qfunc ( sqrt(2*10.^((Eb_No)/10)) )), 'r--' );  % uncoded 2PAM graph
xlabel('Eb/No [dB]'), ylabel('BER')
hold on

semilogy( Eb_No, 1- (1-qfunc ( sqrt(2*10.^(Eb_No/10)))).^4, 'b--' );      % uncoded 2PAM FER graph
xlabel('Eb/No [dB]'), ylabel('BER, FER');
hold on

scatter(Eb_of_No_dB, FER, 20, 'magenta', "X");
hold on
scatter(Eb_of_No_dB, BER, 20, 'green', 'filled');
grid,axis([-2 14 10^(-7) 1]);
legend('Uncoded BER result', 'Uncoded FER result', 'Simulation FER', 'Simulation BER');





