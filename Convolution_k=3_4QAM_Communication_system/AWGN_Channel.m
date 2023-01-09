% AWGN channel

function received_signal = AWGN_channel(code, sigma_v)         % AWGN channel function
    received_signal = code + randn(1, length(code))* sigma_v;    % Y = X + N
end