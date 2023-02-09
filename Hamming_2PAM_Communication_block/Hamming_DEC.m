% Hamming decoding code
% Syndrome decoding

function result = Hamming_DEC(c)
    H = [0 0 0 1 1 1 1;                         % parity-check matrix of (7,4) Hamming code
            0 1 1 0 0 1 1;
            1 0 1 0 1 0 1];

    parity_chk = mod(H*c, 2);              % mod operation bc it's binary
    
    if parity_chk == [0 0 0]                  % no parity error 
        result = [c(3), c(5), c(6), c(7)];
        %disp('no error');
    
    else
        syndrome = 2^0 * parity_chk(3) + 2^1 * parity_chk(2) + 2^2 * parity_chk(1);  % if parity error exist
        c(syndrome) = mod(c(syndrome) + 1, 2);
        result = [c(3), c(5), c(6), c(7)];                                                                          % message bits
        %error = ['e' num2str(syndrome), ' position error occured'];
    end
end