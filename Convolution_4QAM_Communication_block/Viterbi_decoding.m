% Viterbi decoding
% 여기서는 예제로 96bit(=12B)의 메시지를 보낸다고 가정.
function decoded_output = Viterbi_decoding(demodulated_output, num_message_bit)  
    decoded_output = zeros(1, 102);  % 96bits + 6bit(tail bits)
%---------------------------------------------------------------------------------
    % branch metric 계산을 위해서 각 state가 0과 1이 입력으로 들어왔을 때 어떤 출력을 하는지 미리 계산
    binVal = de2bi(0:63, 6, 'left-msb');  
    output_zero = zeros(64,2);
    output_one = zeros(64,2);
    for i = 1:64
        output_zero(i,:) = Convolution_code_not_tail(0, binVal(i,1:6));
        output_one(i,:) = Convolution_code_not_tail(1, binVal(i,1:6));
    end
%--------------------------------------------------------------------------------------
    Path_metric = inf(64, num_message_bit+7);        % 64 X 103 Path metric 저장 위한 matrix
    Path_metric(1, 1) = 0;                                          % 맨 처음 state 
    Message_bit = inf(64, num_message_bit + 7);           % decoding된 메시지를 담을 array
    for i = 1 : num_message_bit+6                           % tail bits 존재
        codeword = demodulated_output(1, 2*i-1: 2*i);  % input을 두개씩 끊어서 codeword 구성
        if i == 1
            Prev_Path_metric_from_zero = zeros(64,1);
            Prev_Path_metric_from_one = zeros(64,1);
            disp("hi")
        else
            for j = 1 : 64
                Current_path = binVal(j,:);                        % ex) current_path = [0 0 1 1 0 0]
                Prev_Path_metric_from_zero = polyval([Current_path(2:6), 0], 2)+1;  % prev_path = [0 1 1 0 0 0]
                Prev_Path_metric_from_one  = polyval([Current_path(2:6), 1], 2)+1;  % prev_path = [0 1 1 0 0 1]

                % BM = min{ 0에서 온 PM + BM, 1에서 온 PM + BM}
                BM_zero = polyval( mod(output_zero(Prev_Path_metric_from_zero) - codeword, 2), 2);
                BM_one  = polyval( mod(output_one(Prev_Path_metric_from_one) - codeword, 2), 2);
                [Path_metric(j, i), Message_bit(j,i)] = min([Prev_Path_metric_from_zero + BM_zero, ...
                    Prev_Path_metric_from_one + BM_one]);
            end
        end
    end
    Message_bit = Message_bit - 1;   % 어떤 인풋을 받아서 현재 state가 되었는지 알기 위한 matrix. Index는 양수라 -1해서 0인지 1인지 저장.
    
%----------------------------------------------------------------------------------
% back tracing
    for i = 1:102
        [~, I] = min(Path_metric(:,i));
        decoded_output(103-i) = Message_bit(I);
    end
end