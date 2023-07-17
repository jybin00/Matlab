% AWGN channel

function received_signal = AWGN_Channel(code, sigma_n)         % AWGN channel function
    received_signal = code + sigma_n*(randn(length(code),1));    % Y = X + N
end