clc
clear
tic
N_m_bit = 30;                                      % Number of message bit
N_frame = 3e5;                                  % Number of frame
test_bit = randsrc(N_frame, N_m_bit, [0 1]);       % test bit generation

Eb_No_dB = (0: 1: 9)';
snr_db = Eb_No_dB;

FER = zeros(1, length(Eb_No_dB));
BER = zeros(1, length(Eb_No_dB));
error = zeros(1, length(Eb_No_dB));
N_f_sim(1,1:length(Eb_No_dB)) = N_frame;

%mode
mode = 2;

trellis = poly2trellis(6, [65, 57]);
%trellis = poly2trellis(3, [7, 5]);
tail_bit = repelem(0, log2(trellis.numStates));
n_mem = trellis.numInputSymbols;

% trellis 정보를 이용해서 입력에 대한 output 구하기.
outputs = trellis.outputs;
output_zero = int2bit(outputs(:,1)', trellis.numInputSymbols)';
output_one  = int2bit(outputs(:,2)', trellis.numInputSymbols)';

for i = 1: length(Eb_No_dB)
    fprintf("Eb/No = %d dB ", Eb_No_dB(i));

    sigma = sqrt(10^(-Eb_No_dB(i)/10));
    for j = 1 : N_frame                             
        % encoded input = (message_bit + tail bit) * 2
        input = test_bit(j, :);                  % tail bit x (pure information bit)
        encoded_input = convenc([input, tail_bit], trellis); % tail bit o
        encoded_input = encoded_input';
        modulated_output = 2 * encoded_input - 1;      % 2PAM or BPSK

        % Signal transmitt through AWGN channel with noise variance sigma_v
        received_signal = awgn(modulated_output, snr_db(i));

        % demodulation
        demodulated_output = received_signal > 0;

        if mode == 1
            decoding = Viterbi_decoding(demodulated_output, N_m_bit);
        elseif mode == 2
            decoding = Vit_gen_hard_dec_mex(demodulated_output, trellis, output_zero, output_one);
        end
        a = nnz(input-decoding);
        if a > 0
            FER(1, i)= FER(1, i) + 1;
            error(1, i) = error(1, i) + a;
        end
        if error(1, i) > 2000
            N_f_sim(1, i) = j;
            break
        end
    end
    BER(1, i) = error(1, i) / (N_f_sim(1, i) * N_m_bit);
    fprintf("BER = %f ", BER(1,i));
    FER(1, i) = FER(1, i) / N_f_sim(1, i);
    fprintf("FER = %f \n", FER(1,i));
end
% FER = FER / N_frame;
toc
%% BER graph
close all
figure(1)
hold on
grid on

axis([0 12 10^-5 1]) 
xticks(0:2:12)

Eb_of_No_db = 0:0.1:12;

linewidth = 1.5;
fontsize = 14;
markersize = 10;

set(gca, "FontName", "Helvativa", "FontSize", fontsize);
set(gca, "yscale", "log");

xlabel("Eb/No", ...
    "FontWeight", "bold");
ylabel("BER", ...
    "FontWeight", "bold");

Eb_of_No_db = -1:0.1:12;

% theorical BER
p_uncoded = semilogy(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), ...
    "color", "#ff0000", ...
    'LineWidth', 1, ...
    'LineStyle', '--');

% sim BER
p_hard = plot(Eb_No_dB, BER, ...
    'Color', '#000000', ...
    'linewidth', 1, ...
    'linestyle', '-', ...
    'Marker', '+', ...
    'MarkerSize', markersize);

lgd = legend([p_uncoded, p_hard], ...
    {'Uncoded 2PAM BER', 'Hard Viterbi v = 5, m = 30'});
lgd.FontSize = fontsize;
lgd.Location = 'best';

%% FER graph

figure(2)
hold on
grid on
axis([0 14 10^-5 10])
xticks(0:2:14)

linewidth = 1.5;
fontsize = 14;
markersize = 10;

Eb_of_No_db = 0:0.1:12;

set(gca, 'FontName', 'Helvatica', 'FontSize', fontsize);
set(gca, 'yscale', 'log');

xlabel("Eb/No", ...
    'FontWeight', 'bold'); 
ylabel('FER', ...
    'FontWeight', 'bold');

% theorical FER
p_uncoded = plot(Eb_of_No_db, 1- (1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) )).^(N_m_bit), ...
    'color', '#ff0000', ...
    'LineWidth', 1, ...
    'LineStyle', '--');

% sim FER
p_hard = plot(Eb_No_dB, FER, ...
    'color', '#000000', ...
    'LineWidth', 1, ...
    'Linestyle', '-', ...
    'marker', 'x', ...
    'markersize', markersize);


lgd = legend([p_uncoded, p_hard], ...
    {'Uncoded 2PAM FER', 'Hard Viterbi FER v = 2, m = 25'});