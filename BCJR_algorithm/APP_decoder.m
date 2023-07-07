% R = 1/2일 때 적용되는 code

function decoded_bit = APP_decoder(input_bit, n_frame, constraintLength, CodeGenerator)

    n_state = 2^(constraintLength-1);

    generator_1 = dec2bin(CodeGenerator(1));
    generator_1 = dec2bin(CodeGenerator(2));

    LLR = zeors(1, length(n_frame));

    alpha = 
    beta = 
    gamma = 

    LLR(i) = log(alpha*gamma*beta / alpha*gamma*beta)

end