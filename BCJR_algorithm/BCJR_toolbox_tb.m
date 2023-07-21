% BCJR test bench
% R = 1/2
tic
clc
close
clear

%number of frame
n_info_bit = 25;
n_frame = 4e5;

%Eb/No db range
Eb_No_dB = (0:10);
SNR_dB = Eb_No_dB;

%trellis
trellis = poly2trellis(3, [6 7]);

%BER, FER variable
FER = zeros(1, length(Eb_No_dB));
BER = zeros(1, length(Eb_No_dB));

%Decoder
BCJR_Decoder = comm.APPDecoder(...
    'TrellisStructure', poly2trellis(3, [6 7]), ...
    'Algorithm', 'True APP', ...
    'CodedBitLLROutputPort', false);

for i = 1:length(Eb_No_dB)
    disp(i)
    for j = 1:n_frame
        info_bit = randi([0 1], 1, n_info_bit);
        messages = convenc(info_bit, trellis);
        % BPSK modulation
        modulated_bit = 2*messages -1;
        % AWGN channel 
        received_bit = awgn(modulated_bit, SNR_dB(i));
        % demodulation
        % demodulated_bit = received_bit > 0;
    
        decoded_bit = BCJR_Decoder(...
            zeros(length(info_bit), 1), (2*received_bit*10^(SNR_dB(i)/10))');
        harddecoded_bit = decoded_bit' > 0;
        BER(i) = BER(i) + nnz(harddecoded_bit - info_bit);
        frame_error = double(nnz(harddecoded_bit ~= info_bit) > 0);
        FER(i) = FER(i) + frame_error;
    end
end

BER = BER/(n_info_bit * n_frame);
FER = FER/n_frame;

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
    "Uncoded 2PAM BER", "APP decoder BER m = 25");

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
    'linewidth', 1, ...
    'linestyle', '--' );

% sim FER
p_hard = plot(Eb_No_dB, FER, ...
    'color', '#000000', ...
    'linewidth', 1, ...
    'linestyle', '-', ...
    'marker', 'x', ...
    'markersize', markersize);


lgd = legend([p_uncoded, p_hard],...
    "Uncoded 2PAM FER", "APP decoder FER m = 25");
lgd.FontSize = fontsize;
lgd.Location = 'best';



