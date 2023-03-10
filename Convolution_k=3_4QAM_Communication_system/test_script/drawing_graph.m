% BER graph
clc
close all
clear

Eb_of_No_db = -1:0.1:15;

F1 = openfig('m = 50/Hard BER.fig');
L1 = findobj(gca,'Type','line');
X1 = L1(1).XData;
Y1 = L1(1).YData;
close

F2 = openfig('m = 25/Hard Viterbi v = 2, m = 25 BER.fig');
L2 = findobj(gca,'Type','line');
X2 = L2(1).XData;
Y2 = L2(1).YData;
close

% theorical BER
figure(1), semilogy(Eb_of_No_db, qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) ), 'r--' );
hold on

plot(X1, Y1, 'kx-')
hold on

plot(X2, Y2, 'bo-')

axis([0 14 0.5*10^-6 1])
xticks(0:2:14)
grid on

xlabel("Eb/No"); 
ylabel('BER');

legend('Uncoded 4QAM BER', 'Hard Viterbi v = 2, m = 50', 'Hard Viterbi v = 2, m = 25')

%%
% FER graph

clc
close all
clear

Eb_of_No_db = -1:0.1:15;

F1 = openfig('m = 50/Soft FER.fig');
L1 = findobj(gca,'Type','line');
X1 = L1(1).XData;
Y1 = L1(1).YData;
close

F2 = openfig('m = 25/Soft Viterbi v = 2, m = 25 FER.fig');
L2 = findobj(gca,'Type','line');
X2 = L2(1).XData;
Y2 = L2(1).YData;
close

% theorical FER
semilogy(Eb_of_No_db, 1- (1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) )).^(50), 'g--' );
hold on
semilogy(Eb_of_No_db, 1- (1-qfunc(sqrt(2*10.^(Eb_of_No_db/10) ) )).^(25), 'r--' );
hold on

plot(X1, Y1, 'kx-')
hold on

plot(X2, Y2, 'bo-')

axis([0 14 10^-5 10])
xticks(0:2:14)
grid on

xlabel("Eb/No"); 
ylabel('FER');

legend('Uncoded 4QAM FER m=50','Uncoded 4QAM FER m=25', 'Soft Viterbi v = 2, m = 50', 'Soft Viterbi v = 2, m = 25')