% Hamming encoder code (7,4,3)_2

function generated_code = Hamming_Enc(input)
    G = [1 1 1 0 0 0 0;
            1 0 0 1 1 0 0;
            0 1 0 1 0 1 0;
            1 1 0 1 0 0 1]; % generating matrix
    generated_code = mod(G'*input, 2);  % result of genearting matrix
end

% generated_code has original codes (c3, c5, c6, c7) = (u1, u2, u3, u4)
% (c1, c2, c4) = (u1+u2+u4, u1+u3+u4, u2+u3+u4)




