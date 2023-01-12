% AWGN channel

function received_signal = AWGN_Channel(code, sigma_v)         % AWGN channel function
    received_signal = code + sigma_v*(randn(1, length(code)) +1i*randn(1, length(code)));    % Y = X + N
end