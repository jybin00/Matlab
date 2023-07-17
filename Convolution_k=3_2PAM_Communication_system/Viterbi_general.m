% Viterbi decoding
% Convolution code는 다음과 같음. (3, [7, 5]) = (7, [111], [101])
function decoded_output = Viterbi_general(demodulated_output, trel)  
    
    tail_bit = trel.numInputSymbols;
    n_mem = trel.numInputSymbols;
    n_states = trel.numStates;
    num_message_bit = length(demodulated_output)/2;
    decoded_output = zeros(1, num_message_bit + tail_bit);  % demodulated_bit + (tail bits)
%---------------------------------------------------------------------------------
    % branch metric 계산을 위해서 각 state가 0과 1이 입력으로 들어왔을 때 어떤 출력을 하는지 미리 계산
    output_zero =  [0 0;  1 1;  1 0;  0 1];

    output_one =  [1 1;  0 0;  0 1;  1 0];  % 1이 입력으로 들어왔을 때 p1, p2
    codeword = reshape(demodulated_output, [n_mem, num_message_bit + tail_bit]);       % input을 두개씩 끊어서 codeword 구성
    codeword = codeword';
%--------------------------------------------------------------------------------------
    Path_metric = inf(n_states, num_message_bit + tail_bit + 1);           % 4 X num_message_bit + 2 + 1 Path metric 저장 위한 matrix
    Message_bit = inf(n_states, num_message_bit + tail_bit + 1);           % path message array
    Survivor_path = zeros(n_states, num_message_bit + tail_bit + 1); 
    for t = 1 : num_message_bit + tail_bit + 1                                      % tail bits 존재
        if t == 1
            Path_metric(1, 1) = 0;	                                                % 맨 처음 state
        elseif t > 1 && t <= tail_bit                                               % t = 2는 state가 펼쳐지는 시간 
            for j = 1 : 2^(t-1)
                Current_state = (n_states/2^(t-1))*(j-1);                     % ex) current_path   = [0 0], [1 0]
                Prev_state_zero = mod(2*Current_state, n_states)+1;            % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(2*Current_state + 1, n_states) +1;       % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_state < 2                                                  % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = nnz(output_zero(Prev_state_zero, :) - codeword(t-1,:));  % 끝자리가 0인 코드에 0이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = nnz(output_zero(Prev_state_one,  :) - codeword(t-1, :));  % 끝자리가 1인 코드에 0이 들어왔을 때
                    [Path_metric(1+Current_state, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                                                     Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(1 + Current_state, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end              
                else	                                                              % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = nnz(output_one(Prev_state_zero, :) - codeword(t-1, :));
                    BM_one  = nnz(output_one(Prev_state_one, :) - codeword(t-1, :));
                    [Path_metric(1+Current_state, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                                                     Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(1+Current_state, t) = 1;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end              
                end
            end
        elseif t > num_message_bit + 1	       % 마지막 2 state는 접히는 시간
            for j = 1 : 2^(num_message_bit + trel.numInputSymbols + 1 -t)			
                Current_state = j-1;                                               % ex) current_path   = [0 0]
                Prev_state_zero = mod(2*Current_state, n_states) + 1;        % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(2*Current_state+1, n_states) + 1;      % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_state < trel.numInputSymbols                            % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = nnz(output_zero(Prev_state_zero, :) - codeword(t-1, :));
                    BM_one  = nnz(output_zero(Prev_state_one,  :) - codeword(t-1, :));
                    [Path_metric(j, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                             Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(j, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end                               
                else                                                                    % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = nnz(output_one(Prev_state_zero, :) - codeword(t-1, :));
                    BM_one  = nnz(output_one(Prev_state_one,  :) - codeword(t-1, :));
                    [Path_metric(j, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                             Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(j, t) = 1;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end                    
                end
            end
        else
            for j = 1 : n_states
                Current_state = j-1;                                 % ex) current_path   = [0 0]
                Prev_state_zero = mod(2*Current_state, n_states)+1;          % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(2*Current_state+1, n_states)+1;        % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_state < trel.numInputSymbols                            % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = nnz(output_zero(Prev_state_zero, :) - codeword(t-1, :));
                    BM_one  = nnz(output_zero(Prev_state_one,  :) - codeword(t-1, :));
                    [Path_metric(j, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                       Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(j, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end

                else                            % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = nnz(output_one(Prev_state_zero, :) - codeword(t-1, :));
                    BM_one  = nnz(output_one(Prev_state_one,  :) - codeword(t-1, :));
                    [Path_metric(j, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                        Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(j, t) = 1;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end    
                end
            end
        end
    end
%----------------------------------------------------------------------------------
% back tracing
    current_state = 00 + 1;                                              % 맨 마지막 state = [0 0]
    %fprintf('%d  ',current_state) 
    Survivor_path = Survivor_path(:, trel.numInputSymbols:end);
    Message_bit = Message_bit(:, trel.numInputSymbols:end);
    for t = num_message_bit + tail_bit :-1 : 1
        decoded_output(1,t) = Message_bit(current_state, t);
        current_state = Survivor_path(current_state, t);
        %fprintf('%d  ',current_state)
    end
    %disp('\n')
    %disp(fliplr(Path_metric));
    %disp(fliplr(Survivor_path));
    decoded_output = decoded_output(1:num_message_bit);
end
