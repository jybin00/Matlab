% BCJR test bench
tic
clc
close
clear

%number of frame and bit
n_info_bit = 25;
n_frame = 1e4;

%Eb/No db range
Eb_No_dB = (0:10);
SNR_dB = Eb_No_dB;
sigma = 10.^(-SNR_dB./10);

%trellis
trellis = poly2trellis(3, [6 7]);

%mode 2 == MAT exchange
%mode 1 == APP decoder
test_mode = 2;

%Decoder
BCJR_Decoder = comm.APPDecoder(...
    'TrellisStructure', poly2trellis(3, [6 7]), ...
    'Algorithm', 'True APP', ...
    'CodedBitLLROutputPort', false);

%BER, FER variable
FER = zeros(1, length(Eb_No_dB));
BER = zeros(1, length(Eb_No_dB));

for i = 1:length(Eb_No_dB)
    disp(i)
    decoded_bit = zeros(1, length(n_info_bit));
    for j = 1:n_frame

        % information bit generation
        info_bit = randi([0 1], 1, n_info_bit);

        %message bit generation
        messages = convenc(info_bit, trellis);

        % BPSK modulation
        modulated_bit = 2*messages-1;
        % AWGN channel 
        received_bit = awgn(modulated_bit, SNR_dB(i));
    
        if test_mode == 1
            decoded_bit = BCJR_Decoder(...
                zeros(length(info_bit), 1),  (2*received_bit*10^(SNR_dB(i)/10))');
    
        elseif test_mode == 2
            decoded_bit = BCJR_matexchange(...
                received_bit, trellis, 10^(-SNR_dB(i)/10) );
        end

        decoded_bit = decoded_bit > 0;
        decoded_bit = 1*reshape(decoded_bit, 1, length(decoded_bit));
        BER(i) = BER(i) + nnz(decoded_bit - info_bit);
    end
    BER(i) = BER(i)/ (n_info_bit * n_frame);
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
lgd.FontSize = fontsize;
lgd.Location = 'best';


