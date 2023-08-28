 

% BCJR test bench
tic
% clc
% close
% clear

%number of frame and bit
n_info_bit = 25;
n_frame = 3e6;
info_bit = randsrc(n_frame, n_info_bit, [0 1]);


%trellis infomation
constraint_length = 3;
generator_polynomial = [5 7];
trellis = poly2trellis(constraint_length, generator_polynomial);

%Eb/No db range & sigma
rate = 1/(constraint_length -1);
Eb_No_dB = (0:7);
SNR_dB = Eb_No_dB + 10*log10(2*rate);
sigma_sq = 10.^(-SNR_dB./10);

trellis_str = string(constraint_length) +' ' +'[' + strjoin(string(generator_polynomial)) + ']' ;

tail_bit = repelem(0, log2(trellis.numStates));
n_mem = log2(trellis.numStates);

% decoding type
% 'True BCJR', 'log-max-MAP', 'log-MAP'
% 'APP BCJR',  'APP max-log', 'APP log'
decoding_mode = "True BCJR";

%BER, FER variable
FER = zeros(1, length(Eb_No_dB));
BER = zeros(1, length(Eb_No_dB));
bit_error = zeros(1, length(Eb_No_dB));
frame_error = zeros(1, length(Eb_No_dB));

for i = 1:length(Eb_No_dB)
    
    % decoder 선언
    BCJR_decoder = F_BCJR_general(trellis, 10^(SNR_dB(i)/10), decoding_mode);

    % simulation 정보
    fprintf("Eb/No = %d dB testing... \n", Eb_No_dB(i));

    decoded_bit = zeros(1, length(n_info_bit));
    tmp_n_frame = n_frame;

    for blk = 1:n_frame

        % Convolution encoding
        messages = convenc([info_bit(blk,:) tail_bit], trellis);

        % BPSK modulation
        modulated_bit = 2*messages-1;
        % AWGN channel 
        received_bit = modulated_bit + sqrt(sigma_sq(i)) * randn(1, length(modulated_bit));
    
        % BCJR decoding
        decoded_bit = BCJR_decoder(zeros(n_info_bit + n_mem, 1), 2*10^(SNR_dB(i)/10)*received_bit');
    
        decoded_bit = 1*reshape(decoded_bit > 0, 1, length(decoded_bit));

        % tail bit 자르기
        decoded_bit = decoded_bit(1:n_info_bit);

        % BER 추가
        bit_error(i) = bit_error(i) + nnz(decoded_bit - info_bit(blk,:));
        f_error = double(nnz(decoded_bit ~= info_bit(blk,:)) > 0);
        frame_error(i) = frame_error(i) + f_error;

        if mod(blk,1000) == 0
            BER(i) = bit_error(i) / (n_info_bit * blk);
            FER(i) = frame_error(i) / blk;
            fprintf(" Eb/No = %.1f [dB]  block = %d/%.2e \n", Eb_No_dB(i), blk, n_frame)
            fprintf(" error : %d\n", bit_error(i))
            fprintf(" BER : %.3e FER : %.3e \n\n", BER(i), FER(i));
        end

        if bit_error(i) > 1000
            BER(i) = bit_error(i) / (n_info_bit * blk);
            FER(i) = frame_error(i) / blk;
            fprintf(" Eb/No = %.1f [dB]  block : %.4e / %.2e \n", Eb_No_dB(i), blk, n_frame)
            fprintf(" error : %d\n", bit_error(i))
            fprintf(" BER: %.3e   FER: %.3e \n\n", BER(i), FER(i));
            break
        end
    end
end

disp("test is done.")

toc
%% BER plot
close
figure(1)
hold on
grid on

axis([0 12 10^-5 1])
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
p_hard = plot(Eb_No_dB, BER, ...
    'color', '#000000', ...
    'linewidth', linewidth, ...
    'linestyle', '-', ...
    'marker', 'x', ...
    'markersize', markersize);



lgd = legend([p_uncoded, p_hard],...
    "Uncoded 2PAM", string(decoding_mode) + " BER m = " + string(n_info_bit));
lgd.FontSize = fontsize;
lgd.Location = 'best';

trellis_info = annotation('textbox', ...
    [.75 .90 .1 .1], ...
    'String',trellis_str,'EdgeColor','red', ...
    'FitBoxToText', 'on');


%% FER plot
figure(2)
hold on
grid on

axis([0 12 10^-4 1])
xticks(0:2:12)

linewidth = 1.5;
fontsize = 14;
markersize = 10;

set(gca, "FontName", "Helvatica", "FontSize", fontsize)
set(gca, "yscale", "log");

xlabel("Eb/No", ...
    "FontWeight", 'bold');
ylabel("FER", ...
    "FontWeight", 'bold');

Eb_of_No_db = -1:0.1:12;

% theorical uncoded FER
p_uncoded = plot(Eb_of_No_db, 1-(1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) )) ).^n_info_bit , ...
    'color', '#ff0000', ...
    'linewidth', linewidth, ...
    'linestyle', '--' );

% sim FER
p_hard = plot(Eb_No_dB, FER, ...
    'color', '#000000', ...
    'linewidth', linewidth, ...
    'linestyle', '-', ...
    'marker', 'x', ...
    'markersize', markersize);

lgd = legend([p_uncoded, p_hard],...
    "Uncoded 2PAM FER", string(decoding_mode) + " FER m = " + string(n_info_bit));
lgd.FontSize = fontsize;
lgd.Location = 'best';

trellis_info = annotation('textbox', ...
    [.75 .90 .1 .1], ...
    'String',trellis_str,'EdgeColor','red', ...
    'FitBoxToText', 'on');

% 외장 모니터 기준
movegui(gca, [1305 530])

