% repetition code test bench
clc
clear
close

eb_No_db = (0:14)';
rate = 1/3;               % 3 repetition
l = 1;                        % bpsk 
snr_db = eb_No_db + 10*log10(rate*l);

N_frame = 100000;
message = randi([0 1], N_frame,1);
bit_error = zeros(1,length(eb_No_db));

for i = 1:length(eb_No_db)
    No = 10^(-snr_db(i)/10)
    code = 2*repmat(message, 1, 3) -1;   % bpsk modulation
    for j = 1:N_frame
        code_k = code(j,:);
        received_signal = awgn(code_k, snr_db(i), 'measured');
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
%%

grid on
semilogy(eb_No_db, ber, '-bo')
hold on
semilogy(eb_No_db, berawgn(eb_No_db, "psk", 2, 'diff'), 'r')
axis([0 14 10^-5 1])
legend('repetition code ber', 'bpsk ber')

xlabel('ber')
ylabel('Eb/No [dB]')
