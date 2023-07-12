% R = 1/2일 때 적용되는 code

function decoded_bit = BCJR_decoder(input_bit, n_frame, constraintLength, SNR)

    n_state = 2^(constraintLength-1);
    n_stage = n_state + 1;

    generator_1 = dec2bin(CodeGenerator(1));
    generator_1 = dec2bin(CodeGenerator(2));

    LLR = zeors(1, length(n_frame));

    % stage는 하나 더 있어야 함. 예를 들어 입력이 2개면 stage는 3개
    alpha = zeors(n_state, n_stage);
    beta = zeros(n_state, n_stage);
    gamma = zeros(n_state, 2* n_stage, constraintLength - 1);

    % gamma 계산을 위해서 각 state가 0과 1이 입력으로 들어왔을 때 어떤 출력을 하는지 미리 계산
    output_zero =  [-1 -1;    1   1;    1 -1;  -1   1];
    output_one =   [  1   1;  -1 -1;  -1   1;    1 -1];

    %input x (두 개씩 끊어서 배열하기)
    input_x = reshape(input_bit, 2, n_frame);

    % 초기 값 설정
    alpha_i = 1;
    beta(1, n_frame) = 1;

    for i = 1:n_frame
        gamma(0+1, i, 1) = exp(4*(SNR/2)*input_x(i) .* output_zero(1)); % (0, 0)
        gamma(1+1, i, 1) = exp(4*(SNR/2)*input_x(i) .* output_zero(2)); % (1, 0)
        gamma(2+1, i, 1) = exp(4*(SNR/2)*input_x(i) .* output_zero(3)); % (2, 1)
        gamma(3+1, i, 1) = exp(4*(SNR/2)*input_x(i) .* output_zero(4)); % (3, 1)

        gamma(0+1, i, 2) = exp(4*(SNR/2)*input_x(i) .* output_one(1)); % (0, 2)
        gamma(1+1, i, 2) = exp(4*(SNR/2)*input_x(i) .* output_one(2)); % (1, 2)
        gamma(2+1, i, 2) = exp(4*(SNR/2)*input_x(i) .* output_one(3)); % (2, 3)  
        gamma(3+1, i, 2) = exp(4*(SNR/2)*input_x(i) .* output_one(4)); % (3, 3)
    end

    % initial stage
    for j = 1 : n_state/2
        alpha_sum = alpha_i * gamma(0+1, i, 1) + alpha_i * gamma(0+1, i, 2);
        alpha(0+1, i) = gamma(0+1, i, 1) / alpha_sum;
        alpha(0+1, )

    LLR(i) = log(alpha(i-1)*gamma(i)*beta(i) / alpha(i-1)*gamma(i)*beta(i));

    end
end





