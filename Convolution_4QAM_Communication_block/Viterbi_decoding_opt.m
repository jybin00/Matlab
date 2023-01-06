% Viterbi decoding
% 여기서는 예제로 96bit(=12B)의 메시지를 보낸다고 가정.
function decoded_output = Viterbi_decoding_opt(demodulated_output, num_message_bit)  
    decoded_output = zeros(1, 102);  % 96bits + 6bit(tail bits)
%---------------------------------------------------------------------------------
    % branch metric 계산을 위해서 각 state가 0과 1이 입력으로 들어왔을 때 어떤 출력을 하는지 미리 계산
    output_zero =     [0 0;  1 1;  0 1;  1 0;  0 0;  1 1;  0 1;  1 0;
                                1 1;  0 0;  1 0;  0 1;  1 1;  0 0;  1 0;  0 1; 
                                1 1;  0 0;  1 0;  0 1;  1 1;  0 0;  1 0;  0 1;
                            	0 0;  1 1;  0 1;  1 0;  0 0;  1 1;  0 1;  1 0;
                                1 0;  0 1;  1 1;  0 0;  1 0;  0 1;  1 1;  0 0;
                            	0 1;  1 0;  0 0;  1 1;  0 1;  1 0;  0 0;  1 1;	
                                0 1;  1 0;  0 0;  1 1;  0 1;  1 0;  0 0;  1 1;
                            	1 0;  0 1;  1 1;  0 0;  1 0;  0 1;  1 1;  0 0];
  
    output_one = [ 1 1;  0 0;	1 0;	0 1;	1 1;	0 0;	1 0;	0 1;
                            0 0;  1 1;	0 1;	1 0;	0 0;	1 1;	0 1;	1 0;
    	                    0 0;  1 1;	0 1;	1 0;	0 0;	1 1;	0 1;	1 0;
                            1 1;  0 0;	1 0;	0 1;	1 1;	0 0;	1 0;	0 1;
                            0 1;  1 0;	0 0;	1 1;	0 1;	1 0;	0 0;	1 1;
                            1 0;  0 1;	1 1;	0 0;	1 0;	0 1;	1 1;	0 0;
                            1 0;  0 1;	1 1;	0 0;	1 0;	0 1;	1 1;	0 0;
                            0 1;  1 0;	0 0;	1 1;	0 1;	1 0;	0 0;	1 1;];		% 1이 입력으로 들어왔을 때 y1, y2
                        	
%--------------------------------------------------------------------------------------
    Path_metric = inf(64, num_message_bit+7);              % 64 X 103 Path metric 저장 위한 matrix
    Message_bit = inf(64, num_message_bit + 7);           % path message array
    for t = 1 : num_message_bit+6+1                             % tail bits 존재
        if t == 1
            Path_metric(1, 1) = 0;                                        % 맨 처음 state
        elseif t > 1 && t <= 6                                            % t = 2~6 까지는 state가 펼쳐지는 시간
            codeword = demodulated_output(1, 1+2*(t-2): 2*(t-1));    % input을 두개씩 끊어서 codeword 구성
            for j = 1 : 2^(t-1)
                Current_path = (64/2^(t-1))*(j-1);          % ex) current_path   = [0 0 1 1 0 0]
                Prev_Path_zero = mod(2*Current_path, 64)+1;          % prev_path_from 0 = [0 1 1 0 0 0] 끝자리 0
                Prev_Path_one  = mod(2*Current_path + 1, 64) +1;   % prev_path_from 1 = [0 1 1 0 0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_path < 32  % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = nnz( mod(output_zero(Prev_Path_zero, :) - codeword, 2)); % 끝자리가 0인 코드에 0이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = nnz( mod(output_zero(Prev_Path_one, :) - codeword, 2)); % 끝자리가 1인 코드에 0이 들어왔을 때
                    Path_metric((1+(64/2^(t-1))*(j-1)), t) = min(Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one);
                    Message_bit((1+(64/2^(t-1))*(j-1)),t) = 0;
                else % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = nnz( mod(output_one(Prev_Path_zero, :) - codeword, 2));
                    BM_one  = nnz( mod(output_one(Prev_Path_one, :) - codeword, 2));
                    Path_metric((1+(64/2^(t-1))*(j-1)), t) = min(Path_metric(Prev_Path_one, t-1) + BM_one, ...
                        Path_metric(Prev_Path_zero, t-1) + BM_zero);
                    Message_bit((1+(64/2^(t-1))*(j-1)), t) = 1;
                end
            end
        elseif t > 97                                                                           % 마지막 t = 98~103는 state가 접히는 시간
            codeword = demodulated_output(1, 1+2*(t-2): 2*(t-1));    % input을 두개씩 끊어서 codeword 구성
            for j = 1 : 2^(103-t)
                Current_path = j-1;                                                       % ex) current_path   = [0 0 1 1 0 0]
                Prev_Path_zero = mod(2*Current_path, 64)+1;             % prev_path_from 0 = [0 1 1 0 0 0] 끝자리 0
                Prev_Path_one  = mod(2*Current_path + 1, 64) +1;     % prev_path_from 1 = [0 1 1 0 0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_path < 32  % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = nnz( mod(output_zero(Prev_Path_zero, :) - codeword, 2));
                    BM_one  = nnz( mod(output_zero(Prev_Path_one, :) - codeword, 2));
                    Path_metric(j, t) = min(Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one);
                    Message_bit(j,t) = 0;
                else % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = nnz( mod(output_one(Prev_Path_zero, :) - codeword, 2));
                    BM_one  = nnz( mod(output_one(Prev_Path_one, :) - codeword, 2));
                    Path_metric(j, t) = min(Path_metric(Prev_Path_one, t-1) + BM_one, ...
                        Path_metric(Prev_Path_zero, t-1) + BM_zero);
                    Message_bit(j,t) = 1;
                end
            end
        else
            codeword = demodulated_output(1, 1+2*(t-2): 2*(t-1));    % input을 두개씩 끊어서 codeword 구성
            for j = 1 : 64
                Current_path = j-1;                                                         % ex) current_path   = [0 0 1 1 0 0]
                Prev_Path_zero = mod(2*Current_path, 64)+1;               % prev_path_from 0 = [0 1 1 0 0 0] 끝자리 0
                Prev_Path_one  = mod(2*Current_path + 1, 64) +1;        % prev_path_from 1 = [0 1 1 0 0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_path < 32   % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = nnz( mod(output_zero(Prev_Path_zero, :) - codeword, 2));
                    BM_one  = nnz( mod(output_zero(Prev_Path_one, :) - codeword, 2));
                    Path_metric(j, t) = min(Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one);
                    Message_bit(j,t) = 0;
                else % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = nnz( mod(output_one(Prev_Path_zero, :) - codeword, 2));
                    BM_one  = nnz( mod(output_one(Prev_Path_one, :) - codeword, 2));
                    Path_metric(j, t) = min(Path_metric(Prev_Path_zero, t-1) + BM_zero, ...
                        Path_metric(Prev_Path_one, t-1) + BM_one);
                    Message_bit(j,t) = 1;
                end
            end
        end
    end
%----------------------------------------------------------------------------------
% back tracing
    back_path = 0;                                          % 맨 마지막 state = [0 0 0 0 0 0]
    Message_bit = Message_bit(:,2:103);         % Message matrix 한 칸 당기기.
    for t = 1:102
        prev_back_path0 = mod(2*back_path, 64)+1;  
        prev_back_path1 = mod(2*back_path + 1, 64) +1;
        [~, I] = min([Path_metric(prev_back_path0, 103-t), Path_metric(prev_back_path1, 103-t )]);
        if I == 1
            decoded_output(1,103-t) = Message_bit(back_path+1, 103-t);
            back_path = mod(2*back_path, 64);
            % disp(prev_back_path0+1)
        else
            decoded_output(1,103-t) = Message_bit(back_path+1, 103-t);
            back_path = mod(2*back_path + 1, 64);
            % disp(prev_back_path1+1)
        end
    end
    decoded_output = decoded_output(1:96);
end

