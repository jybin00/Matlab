% repetition code test bench
clc
clear
close
tic

%%rng("default")

eb_No_db = (0:15)';
rate = 1/3;               % 3 repetition
l = 1;                        % bpsk 
snr_db = eb_No_db + 10*log10(rate*l);

N_frame = 1000000;
message = randi([0 1], N_frame,1);
bit_error = zeros(1,length(eb_No_db));
No = 10.^(-snr_db/10);
q = zeros(length(No),1);

parfor i = 1:length(eb_No_db)
    received_signal = zeros(1,3);
    sigma = sqrt(No(i)/2);
    code = 2*repmat(message, 1, 3) -1;   % bpsk modulation
    for j = 1:N_frame
        code_k = code(j,:);
        received_signal = code_k + sigma*randn(1,3);
        if sum((received_signal - 1).^2) < sum((received_signal + 1).^2)
            decoding = 1;
            if decoding == code_k(1)
                continue;
            else
                bit_error(1,i) = bit_error(1,i) + 1;
            end
        else
            decoding = -1;
            if decoding == code_k(1)
                continue;
            else
                bit_error(1,i) = bit_error(1,i) + 1;
            end
        end  
    end
end
ber = bit_error/N_frame
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
grid on
semilogy(eb_No_db, 1-q, 'kx-')
hold on
semilogy(eb_No_db, ber, '-bo')
hold on
semilogy(eb_No_db, berawgn(eb_No_db, "psk", 2, 'nondiff'), 'r')
axis([0 13 10^-5 1])
legend('Soft decision theoretical ber', 'Repetition code soft decision ber', 'BPSK ber')

ylabel('BER')
xlabel('Eb/No [dB]')
