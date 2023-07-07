% BCJR test bench
tic
clc
close
clear

%number of frame
n_info_bit = 3e6;

%Eb/No db range
Eb_No_dB = (0:10);
SNR_dB = Eb_No_dB;

%trellis
trellis = poly2trellis(3, [6 7]);
info_bit = randi([0 1], 1, n_info_bit);
messages = convenc(info_bit, trellis);

%BER, FER variable
FER = zeros(1, length(Eb_No_dB));
BER = zeros(1, length(Eb_No_dB));

for i = 1:length(Eb_No_dB)
    % BPSK modulation
    modulated_bit = 2*messages -1;
    % AWGN channel 
    received_bit = awgn(modulated_bit, SNR_dB(i));
    % demodulation
    % demodulated_bit = received_bit > 0;

    decoded_bit = BCJR_Decoder(...
        zeros(length(info_bit), 1), received_bit');
    harddecoded_bit = decoded_bit > 0;
    BER(i) = nnz(harddecoded_bit' - info_bit)/ n_info_bit;
end

disp(BER)

toc
%%
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
    'linewidth', 1, ...
    'linestyle', '--' );

% sim BER
p_hard = plot(Eb_No_dB, BER, ...
    'color', '#000000', ...
    'linewidth', 1, ...
    'linestyle', '-', ...
    'marker', 'x', ...
    'markersize', markersize);


lgd = legend([p_uncoded, p_hard],...
    "Uncoded 2PAM", "hard BCJR");


