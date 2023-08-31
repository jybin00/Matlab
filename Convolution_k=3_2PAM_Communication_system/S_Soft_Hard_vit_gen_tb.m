%% rate = 1/2에 해당하는 모든 convolution code에 대해서 실행되는 vitdec 알고리즘 만들기.
clc
clear
tic
N_m_bit = 1000;                   % Number of message bit
N_frame = 5e6;                   % Number of frame


Eb_No_dB = (0:7)';

% trellis 정보를 이용해서 입력에 대한 output 구하기.
constraint_length = 3;
trellis = poly2trellis(constraint_length + 1, [13, 15, 17]);

output_zero = 2 * int2bit(trellis.outputs(:,1)', log2(trellis.numOutputSymbols))' - 1;
output_one  = 2 * int2bit(trellis.outputs(:,2)', log2(trellis.numOutputSymbols))' - 1;

tail_bit = repelem(0, log2(trellis.numStates));
% test_bit = randsrc(N_frame, N_m_bit, [0 1]);     % test bit generation

rate = 1/constraint_length;
snr_dB = Eb_No_dB + 10*log10(2*rate);
sigma_sq = 10.^(-snr_dB/10);

FER_soft = zeros(1, length(Eb_No_dB));
% FER_hard = zeros(1, length(Eb_No_dB));

BER_soft = zeros(1, length(Eb_No_dB));
% BER_hard = zeros(1, length(Eb_No_dB));

bit_error_soft = zeros(1, length(Eb_No_dB));
frame_error_soft = zeros(1, length(Eb_No_dB));
% error_hard = zeros(1, length(Eb_No_dB));

N_f_sim(1,1:length(Eb_No_dB)) = N_frame;

% hard_decoder = F_Vitdec(trellis, N_m_bit, 'type', 'hard');
soft_decoder = F_Vitdec(trellis, N_m_bit, 'type', 'soft');


for n = 1 : length(Eb_No_dB)         % Eb를 바꿔가면서 계산
    fprintf("Eb of No : %d dB\n", n-1) % 진행상황 확인을 위한 인덱스
    for blk = 1 : N_frame              % frame 개수만큼 계산

        % 1~ message bit || 1+ message bit ~ 2* message bit ....
        % encoded input = (message_bit + tail bit) * 2
        % input = test_bit(blk, :);
        input = randsrc(1, N_m_bit, [0 1]); 
        encoded_input = convenc([input, tail_bit], trellis);  

        % modulation -> channel -> demodulation
        modulated_output = 2*encoded_input -1;

        % Signal transmitt through AWGN channel with noise variance sigma_v
        received_signal = modulated_output + sqrt(sigma_sq(n)) * randn(1, length(modulated_output));

        % demodulation
        demodulated_output = received_signal > 0;
        
        %decoding_soft = soft_decoder(received_signal);
        decoding_soft = F_Vit_gen_soft_dec_mex(received_signal, trellis, output_zero, output_one, 1/3);
        % decoding_hard = hard_decoder(demodulated_output);

        error_s = sum(mod(input+decoding_soft,2));
        % error_h = sum(mod(input+decoding_hard,2));
        if error_s > 0
            frame_error_soft(1, n)= frame_error_soft(1, n) + 1;
            bit_error_soft(1, n) = bit_error_soft(1, n) + error_s;
        end
        % if error_h > 0
        %     FER_hard(1, n)= FER_hard(1, n) + 1;
        %     error_hard(1, n) = error_hard(1, n) + error_h;
        % end
        
        if mod(blk, 3000) == 0
            BER_soft(1, n) = bit_error_soft(1, n) / (blk * N_m_bit);
            % BER_hard(1, n) = error_hard(1, n) / (blk * N_m_bit);
        
            FER_soft(1, n) = frame_error_soft(1, n) / blk;
            % FER_hard(1, n) = FER_hard(1, n) / blk;

            fprintf(" Eb/No = %.1f [dB]  block = %.2e/%.2e \n", Eb_No_dB(n), blk, N_frame);
            fprintf(" error : %d\n", bit_error_soft(1, n))
            fprintf(" BER : %.3e   FER : %.3e \n\n", BER_soft(1, n), FER_soft(1, n));
        end


        if bit_error_soft(1, n) > 10000
            BER_soft(1, n) = bit_error_soft(1, n) / (blk * N_m_bit);
            % BER_hard(1, n) = error_hard(1, n) / (blk * N_m_bit);
        
            FER_soft(1, n) = frame_error_soft(1, n) / blk;
            % FER_hard(1, n) = FER_hard(1, n) / blk;

            fprintf(" Eb/No = %.1f [dB]  block = %.2e/%.2e \n", Eb_No_dB(n), blk, N_frame);
            fprintf(" error : %d\n", bit_error_soft(1, n))
            fprintf(" BER : %.3e   FER : %.3e \n\n", BER_soft(1, n), FER_soft(1, n));
            break
        end
        
    end

end
toc

%% BER graph
close all
figure(1)
hold on
grid on

axis([0 12 10^-6 1])
xticks(0:2:12)

linewidth = 1.5;
fontsize = 14;
markersize = 10;

set(gca, "FontName", "Helvatica", "FontSize", fontsize)
set(gca, "yscale", "log");

xlabel("Eb/No", ...
    "FontWeight", 'bold');
ylabel("BER", ...
    "FontWeight", 'bold');

Eb_of_No_db = -1:0.1:12;

% theorical uncoded BER
p_uncoded = plot(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), ...
    'color', '#ff0000', ...
    'linewidth', linewidth, ...
    'linestyle', '--' );

% sim BER
% p_hard = plot(Eb_No_dB, BER_hard, ...
%     'color', '#000000', ...
%     'linewidth', linewidth, ...
%     'linestyle', '-', ...
%     'marker', 'x', ...
%     'markersize', markersize);
% 
p_soft = plot(Eb_No_dB, BER_soft, ...
    'color', '#0000ff', ...
    'linewidth', linewidth, ...
    'linestyle', '-', ...
    'marker', 'o', ...
    'markersize', markersize);

%legend('Uncoded 2PAM BER', 'Hard wrong Viterbi v = 2, m = 50', 'Hard correct Viterbi v = 2, m = 50')
lgd = legend([p_uncoded, p_soft], ...
    {'Uncoded BPSK BER', 'Soft Viterbi , m = '+string(N_m_bit)});
% lgd = legend([p_uncoded, p_hard, p_soft], ...
%     {'Uncoded BPSK BER', 'Hard Viterbi v = 2, m = 100', 'Soft Viterbi v = 2, m = 100'});
lgd.FontSize = fontsize;
lgd.Location = 'best';

%% FER graph

figure(2)
hold on
grid on
axis([0 12 10^-4 1])
xticks([0:2:12])

linewidth = 1.5;
fontsize = 13;
markersize = 9;

set(gca, 'FontName', "Helvatica", 'FontSize', fontsize)
set(gca, 'yscale', 'log');

xlabel("Eb/No", ...
    'FontWeight', 'bold');
ylabel("FER", ...
    'FontWeight', 'bold');

Eb_of_No_db = -1:0.1:12;

% theorical FER
p_theory = plot(Eb_of_No_db, 1- (1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) )).^(N_m_bit), ...
    'color', '#ff0000', ...
    'linewidth', linewidth, ...
    'linestyle', '--' );

% sim FER
% p_hard = plot(Eb_No_dB, FER_hard, ...
%     'color', '#000000', ...
%     'linewidth', linewidth, ...
%     'linestyle', '-', ...
%     'marker', 'x', ...
%     'markersize', markersize);

p_soft = plot(Eb_No_dB, FER_soft, ...
    'color', '#0000ff', ...
    'linewidth', linewidth, ...
    'linestyle', '-', ...
    'marker', 'o', ...
    'markersize', markersize);


% lgd = legend([p_theory, p_hard, p_soft], ...
%     {'Uncoded 2PAM FER', 'Hard Viterbi v = 2, m = 50', 'Soft Viterbi v = 2, m = 50'});
lgd = legend([p_theory, p_soft], ...
    {'Uncoded BPSK FER', 'Soft Viterbi, m = '+string(N_m_bit)});
lgd.FontSize = fontsize;
lgd.Location = 'best';


% legend('Uncoded 4QAM FER', 'Soft Viterbi v = 2, m = 50')