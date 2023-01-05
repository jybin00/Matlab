% Viterbi decoding
% 여기서는 예제로 96bit(=12B)의 메시지를 보낸다고 가정.
function decoded_output = Viterbi_decoding(demodulated_output, num_message_bit)  
    decoded_output = zeros(1, 102);  % 96bits + 6bit(tail bits)
%---------------------------------------------------------------------------------
    % branch metric 계산을 위해서 각 state가 0과 1이 입력으로 들어왔을 때 어떤 출력을 하는지 미리 계산
    binVal = de2bi(0:63, 6, 'left-msb');  
    output_zero = zeros(64,2);
    output_one = zeros(64,2);
    for t = 1:64
        output_zero(t,:) = Convolution_code_not_tail(0, binVal(t,1:6));
        output_one(t,:) = Convolution_code_not_tail(1, binVal(t,1:6));
    end
%--------------------------------------------------------------------------------------
    Path_metric = inf(64, num_message_bit+7);              % 64 X 103 Path metric 저장 위한 matrix
    Message_bit = inf(64, num_message_bit + 7);           % path message array
    Current_path = zeros(1, 6);
    for t = 1 : num_message_bit+6+1                             % tail bits 존재
        if t == 1
            Path_metric(1, 1) = 0;                                        % 맨 처음 state
        elseif t > 1 && t < 6                                               % t = 2~6 까지는 state가 펼쳐지는 시간
            codeword = demodulated_output(1, 1+2*(t-2): 2*(t-1));    % input을 두개씩 끊어서 codeword 구성
            for j = 1 : 2^(t-1)
                Current_path = binVal(1+(64/2^(t-1))*(j-1), : );          % ex) current_path   = [0 0 1 1 0 0]
                Prev_Path_zero = polyval([Current_path(2:6), 0], 2)+1;  % prev_path_from 0 = [0 1 1 0 0 0] 끝자리 0
                Prev_Path_one  = polyval([Current_path(2:6), 1], 2)+1;  % prev_path_from 1 = [0 1 1 0 0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_path(1) == 0  % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = polyval( mod(output_zero(Prev_Path_zero) - codeword, 2), 2);
                    BM_one  = polyval( mod(output_zero(Prev_Path_one) - codeword, 2), 2);
                    [Path_metric((1+(64/2^(t-1))*(j-1)), t), Message_bit((1+(64/2^(t-1))*(j-1)),t)] = min([Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one]);
                else % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = polyval( mod(output_one(Prev_Path_zero) - codeword, 2), 2);
                    BM_one  = polyval( mod(output_one(Prev_Path_one) - codeword, 2), 2);
                    [Path_metric((1+(64/2^(t-1))*(j-1)), t), Message_bit((1+(64/2^(t-1))*(j-1)),t)] = min([Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one]);
                end
            end
        elseif t > 97                                                                          % 마지막 t = 98~103는 state가 접히는 시간
            codeword = demodulated_output(1, 1+2*(t-2): 2*(t-1));    % input을 두개씩 끊어서 codeword 구성
            for j = 1 : 2^(103-t)
                Current_path = binVal(1+(64/2^(t-1))*(j-1), : );            % ex) current_path   = [0 0 1 1 0 0]
                Prev_Path_zero = polyval([Current_path(2:6), 0], 2)+1;  % prev_path_from 0 = [0 1 1 0 0 0] 끝자리 0
                Prev_Path_one  = polyval([Current_path(2:6), 1], 2)+1;  % prev_path_from 1 = [0 1 1 0 0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_path(1) == 0  % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = polyval( mod(output_zero(Prev_Path_zero) - codeword, 2), 2);
                    BM_one  = polyval( mod(output_zero(Prev_Path_one) - codeword, 2), 2);
                    [Path_metric(1+(64/2^(t-1))*(j-1), t), Message_bit((1+(64/2^(t-1))*(j-1)),t)] = min([Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one]);
                else % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = polyval( mod(output_one(Prev_Path_zero) - codeword, 2), 2);
                    BM_one  = polyval( mod(output_one(Prev_Path_one) - codeword, 2), 2);
                    [Path_metric(1+(64/2^(t-1))*(j-1), t), Message_bit((1+(64/2^(t-1))*(j-1)),t)] = min([Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one]);
                end
            end
        else
            codeword = demodulated_output(1, 1+2*(t-2): 2*(t-1));    % input을 두개씩 끊어서 codeword 구성
            for j = 1 : 64
                Current_path = binVal(j,:);                                             % ex) current_path   = [0 0 1 1 0 0]
                Prev_Path_zero = polyval([Current_path(2:6), 0], 2)+1;  % prev_path_from 0 = [0 1 1 0 0 0] 끝자리 0
                Prev_Path_one  = polyval([Current_path(2:6), 1], 2)+1;  % prev_path_from 1 = [0 1 1 0 0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_path(1) == 0  % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = polyval( mod(output_zero(Prev_Path_zero) - codeword, 2), 2);
                    BM_one  = polyval( mod(output_zero(Prev_Path_one) - codeword, 2), 2);
                    [Path_metric(j, t), Message_bit((1+(64/2^(t-1))*(j-1)),t)] = min([Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one]);
                else % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = polyval( mod(output_one(Prev_Path_zero) - codeword, 2), 2);
                    BM_one  = polyval( mod(output_one(Prev_Path_one) - codeword, 2), 2);
                    [Path_metric(j, t), Message_bit((1+(64/2^(t-1))*(j-1)),t)] = min([Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one]);
                end
            end
        end
    end
    Message_bit = Message_bit - 1;   % 어떤 인풋을 받아서 현재 state가 되었는지 알기 위한 matrix. Index는 양수라 -1해서 0인지 1인지 저장.
    
%----------------------------------------------------------------------------------
% back tracing
    for t = 1:102
        [~, I] = min(Path_metric(:,103-t));
        decoded_output(103-t) = Message_bit(I);
        decoded_output = decoded_output(1:96);
    end
end