function decoder = ML_decoder(codeword_dict, channel)
switch channel
    case 'BSC'
        decoder = @decode_BSC;
    case 'BEC'
        decoder = @decode_BEC;
    case 'AWGN'
        decoder = @decode_AWGN;
    otherwise
        error('Unsupported channel %s', channel);
end
    function x_hat_enum = decode_BSC(y)
        % ML decoding
        [~, idx] = min(sum(mod(y + codeword_dict, 2), 2));
        x_hat_enum = idx - 1;
    end
    function x_hat_enum = decode_BEC(y)
        alived_idx = ~isnan(y);
        codeword_dict_erased = codeword_dict(:, alived_idx);
        [~, idx] = min(sum(mod(y(alived_idx) + codeword_dict_erased, 2), 2));
        x_hat_enum = idx - 1;
    end
    function x_hat_enum = decode_AWGN(y)
        % ML decoding
        [~, idx] = min(sum(abs(y - codeword_dict).^2, 2));
        x_hat_enum = idx - 1;
    end
end