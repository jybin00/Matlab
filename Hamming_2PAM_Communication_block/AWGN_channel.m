% AWGN channel

function received_signal = AWGN_channel(code, sigma_v)         % AWGN channel function
    received_signal = code + randn(length(code), 1)* sigma_v;    % Y = X + N
end