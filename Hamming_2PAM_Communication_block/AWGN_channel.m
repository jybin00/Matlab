% AWGN channel

function received_signal = AWGN_channel(code, noise_v)         % AWGN channel function
    received_signal = code + randn(length(code), 1)* noise_v;    % Y = X + N(=n*sigma)
end