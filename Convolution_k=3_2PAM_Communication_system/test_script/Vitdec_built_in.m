clc
clear
tic
% 매트랩 내장 함수인 vitdec을 이용해서 실제로 제대로 구현한 것인지 확인해보기.

%system parameter
Eb_No_dB = (0:10)';
SNR_dB = Eb_No_dB;
BER_hard = zeros(1, length(Eb_No_dB));
BER_soft = zeros(1, length(Eb_No_dB));
FER_hard = zeros(1, length(Eb_No_dB));
tb = 50;

% input
number_of_msg = 4e6;
msg = randi([0 1], 1, number_of_msg);

% convolution code
trellis = poly2trellis(6,[65 57]);
coded_bit = convenc(msg, trellis);

% mapping
symbol = 2*coded_bit - 1;

for i = 1: length(Eb_No_dB)
    fprintf("Eb/No = %d dB \n", Eb_No_dB(i));
    % channel
    received_bit = awgn(symbol, SNR_dB(i));
    
    % demapping
    demodulated_output = received_bit > 0;
    
    % decoding
    decoding_hard = vitdec(demodulated_output, trellis, tb, 'term', 'hard');
    % llr 계산이 반대로 되어 있어서 - 붙여야 한다ter c/ 
    decoding_soft =  vitdec(-received_bit,     trellis, tb, 'term', 'unquant');
    
    % FER, BER measure
    BER_hard(i) = nnz(decoding_hard - msg) / number_of_msg;
    BER_soft(i) = nnz(decoding_soft - msg) / number_of_msg;
end
toc

%%
close all
figure(1)
hold on
grid on

linewidth = 1.5;
fontsize = 14;
markersize = 10;

set(gca, 'FontName', 'Helvatica', 'FontSize', fontsize)
set(gca, 'yscale', 'log');

xlabel("Eb/No", ...
    'FontWeight', 'bold'); 
ylabel('BER', ...
    'FontWeight', 'bold');


Eb_of_No_db = -1:0.1:12;
% theorical uncoded BER
% semilogy(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), 'r--' );
p_uncoded = plot(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), ...
    'color', '#ff0000', ...
    'linewidth', 1, ...
    'linestyle', '--');

% sim BER
% semilogy(Eb_No_dB, BER_hard, 'kx-')
p_hard_viterbi = plot(Eb_No_dB, BER_hard, ...
    'color', '#000000', ...
    'linewidth', linewidth, ...
    'linestyle', '-', ...
    'marker', 'x', ...
    'markersize', markersize);

% semilogy(Eb_No_dB, BER_soft, 'bo-')
p_soft_viterbi = plot(Eb_No_dB, BER_soft, ...
    'color', '#0000ff', ...
    'linewidth', linewidth, ...
    'linestyle', '-', ...
    'marker', 'o', ...
    'markersize', markersize);

axis([0 12 10^-5 1])
xticks(0:2:12)



%legend('Uncoded 2PAM BER', 'Hard wrong Viterbi v = 2, m = 50', 'Hard correct Viterbi v = 2, m = 50')
lgd = legend([p_uncoded, p_hard_viterbi, p_soft_viterbi], ...
    {'Uncoded 2PAM BER', 'Hard Viterbi v = 5', 'Soft Viterbi v = 5'});
lgd.FontSize = fontsize;
lgd.Location = 'best';


% annotation('doublearrow',[0.448214285714286 0.675],...
%     [0.270428571428571 0.276190476190476]);

% % textbox 생성
% annotation('textbox',...
%     [0.478571428571429 0.204952380952382 0.0959821428571429 0.0654761904761905],...
%     'String',{'5 dB'},...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

