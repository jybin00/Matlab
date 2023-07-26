 

% BCJR test bench
tic
clc
close
clear

%number of frame and bit
n_info_bit = 20;
n_frame = 2e1;
info_bit = randsrc(n_frame, n_info_bit, [0 1]);

%Eb/No db range
Eb_No_dB = (0:5);
SNR_dB = Eb_No_dB;
sigma = 10.^(-SNR_dB./10);

%trellis
constraint_length = [4 4];
generator_polynomial = [17 15 13; 17 15 13];
trellis = poly2trellis(constraint_length, generator_polynomial);
% trellis = poly2trellis([4 4],[15 13 17; 17 15 13]);

trellis_str = string(constraint_length) +' ' +'[' + strjoin(string(generator_polynomial)) + ']' ;

tail_bit = repelem(0, log2(trellis.numStates));
n_mem = log2(trellis.numStates);

% 'True BCJR', 'log-max-MAP', 'log-MAP' 
% 'APP_BCJR',  'APP_max_log', 'APP_log'
decoding_mode = "log-max-MAP";


%BER, FER variable
FER = zeros(1, length(Eb_No_dB));
BER = zeros(1, length(Eb_No_dB));

for i = 1:length(Eb_No_dB)

    BCJR_decoder = BCJR_general(trellis, 10^(SNR_dB(i)/10), decoding_mode);
    fprintf("Eb/No = %d dB testing... \n", i-1);
    decoded_bit = zeros(1, length(n_info_bit));
    tmp_n_frame = n_frame;

    for j = 1:n_frame
        
        messages = convenc([info_bit(j,:) tail_bit], trellis);

        % BPSK modulation
        modulated_bit = 2*messages-1;
        % AWGN channel 
        received_bit = awgn(modulated_bit, SNR_dB(i));
    
        decoded_bit = BCJR_decoder(2*10^(SNR_dB(i)/10)*received_bit);
    
        decoded_bit = 1*reshape(decoded_bit > 0, 1, length(decoded_bit));
        decoded_bit = decoded_bit(1:n_info_bit);
        BER(i) = BER(i) + nnz(decoded_bit - info_bit(j,:));
        frame_error = double(nnz(decoded_bit ~= info_bit(j,:)) > 0);
        FER(i) = FER(i) + frame_error;
        if BER(i) > 4000
            tmp_n_frame = j;
            break
        end
    end
    BER(i) = BER(i) / (n_info_bit * tmp_n_frame);
    FER(i) = FER(i) / tmp_n_frame;
    fprintf(" BER : %.3e \n FER : %.3e \n\n", BER(i), FER(i));
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

txt = text('Position', [8 10^-1], ...
    'String', trellis_str, ...
    'FontSize', 16, ...
    'FontWeight', 'bold');


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

txt = text('Position', [8 10^-1], ...
    'String', trellis_str, ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

% 외장 모니터 기준
movegui(gca, [1305 530])

