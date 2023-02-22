% repetition code test bench
clc
close
tic

%%rng("default")

eb_No_db = (0:15)';
rate = 1/3;               % 3 repetition
k = 1;                        % bpsk 
snr_db = eb_No_db + 10*log10(rate*k);

N_frame = 1000000;
message = randi([0 1], N_frame,1);
bit_s_e = zeros(1,length(eb_No_db));
bit_h_e = zeros(1,length(eb_No_db));
No = 10.^(-snr_db/10);
q = zeros(length(No),1);

parfor i = 1:length(eb_No_db)
    received_signal = zeros(1,3);
    sigma = sqrt(No(i)/2);
    code = 2*repmat(message, 1, 3) -1;   % bpsk modulation
    for j = 1:N_frame
        code_k = code(j,:);
        % received_signal = code_k + sigma*randn(1,3);    % 이렇게 하니까 잘 되는데 왜 awgn로 하면 안되지??
        received_signal = awgn(code_k, 3+snr_db(i))         % 3dB를 더 해줘야 한다. 왜?? 직접 겪는 노이즈는 No/2니까. 
        if sum((received_signal - 1).^2) < sum((received_signal + 1).^2)
            decoding = 1;
        else
            decoding = -1;
        end  

        if decoding ~= code_k(1)
            bit_s_e(1,i) = bit_s_e(1,i) + 1;
        end

        received_signal = received_signal > 0;
        if sum(received_signal) > 1    % numel은 행렬의 원소의 개수를 반환하는 함수.
            decoding = 1;
        else
            decoding = -1;
        end

        if decoding ~= code_k(1)
            bit_h_e(1,i) = bit_h_e(1,i) + 1;
        end
    end
end
ber_s = bit_s_e/N_frame;
ber_h = bit_h_e/N_frame;
toc
%%
close all
tic
parfor i = 1:length(No)
    disp(i)
    fun = @(x,y,z) (1/(pi*No(i)).^(3/2))*exp(-((x-1).^2 + (y-1).^2 + (z-1).^2)/No(i));
    
    xmin = -inf;
    xmax = inf;
    ymin = -inf;
    ymax = inf;
    zmin = @(x,y) -y-x;
    zmax = inf;
    
    q(i,1) = integral3(fun, xmin, xmax, ymin, ymax, zmin, zmax, 'Method', 'iterated');
end
toc
%%
close all

%semilogy(eb_No_db, 1-q, 'mdiamond-')
%hold on
%semilogy(eb_No_db, ber_s, '-bo')
%hold on
semilogy(eb_No_db, ber_h, '-ko')
hold on
semilogy(eb_No_db, 3*qfunc(sqrt((2/3)*10.^(eb_No_db/10))).^2-2*qfunc(sqrt((2/3)*10.^(eb_No_db/10))).^3, 'bx-' )
hold on
semilogy(eb_No_db, berawgn(eb_No_db, "psk", 2, 'nondiff'), 'r')
axis([0 13 10^-5 1])
%legend('Soft decision theoretical ber', '3-Repetition code soft decision ber', 'BPSK ber')
legend( '1/3-rate repetition hard ber', '1/3-drate repetition hard theretical ber','BPSK ber')
grid on

ylabel('BER')
xlabel('Eb/No [dB]')
