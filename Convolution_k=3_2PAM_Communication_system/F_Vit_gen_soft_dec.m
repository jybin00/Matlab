% Viterbi decoding
% Convolution code는 trellis 정보에 따라 달라짐.
% R = 1/2 code
function decoded_output = F_Vit_gen_soft_dec(received_signal, trel, output_zero, output_one, rate)  
    
    tail_bit = log2(trel.numStates);
    n_mem = log2(trel.numStates);
    n_states = trel.numStates;
    n_inputSymbols = trel.numInputSymbols;
    N = 1/rate;
    
    num_message_bit = length(received_signal)/N - tail_bit;
    decoded_output = zeros(1, num_message_bit + tail_bit);  % demodulated_bit + (tail bits)
%---------------------------------------------------------------------------------
    % input을 n_mem개씩 끊어서 received_vector 구성 
    received_vector = reshape(received_signal, [N, num_message_bit + tail_bit]);     
    received_vector = received_vector';
%--------------------------------------------------------------------------------------
    % # of states X (num_message_bit + tail bit + 1) Path metric 저장 위한 matrix
    Path_metric = inf(n_states, num_message_bit + tail_bit + 1);          
    Message_bit = inf(n_states, num_message_bit + tail_bit + 1);           % path message array
    Survivor_path = zeros(n_states, num_message_bit + tail_bit + 1); 
    for t = 1 : num_message_bit + tail_bit + 1  
        % 맨 처음 state
        if t == 1
            Path_metric(1, 1) = 0;	

        % t = n_state까지는 state가 펼쳐지는 시간 
        elseif t > 1 && t <= n_mem                                                 
            for j = 1 : 2^(t-1)
                Current_state = (n_states/2^(t-1))*(j-1);                            % ex) current_path = [0 0], [1 0]
                Prev_state_zero = mod(n_inputSymbols * Current_state, n_states)+1;            % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(n_inputSymbols * Current_state + 1, n_states) +1;       % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}

                % 입력 0을 받아서 현재 값이 되었을 경우
                if Current_state < n_states / n_inputSymbols                         
                    BM_zero = sum((output_zero(Prev_state_zero, :) - received_vector(t-1,:)).^2);  % 끝자리가 0인 코드에 0이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = sum((output_zero(Prev_state_one, :) - received_vector(t-1, :)).^2);  % 끝자리가 1인 코드에 0이 들어왔을 때
                    [Path_metric(1+Current_state, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                                                     Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(1 + Current_state, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end        

                % 입력 1을 받아서 현재 값이 되었을 경우    
                else	                                                            
                    BM_zero = sum((output_one(Prev_state_zero, :) - received_vector(t-1,:)).^2);  % 끝자리가 0인 코드에 1이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = sum((output_one(Prev_state_one, :) - received_vector(t-1, :)).^2);  % 끝자리가 1인 코드에 1이 들어왔을 때
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
        % 마지막은 tail bit 갯수만큼 trellis가 접히는 시간
        elseif t > num_message_bit + 1	       
            for j = 1 : 2^(num_message_bit + n_mem + 1 -t)			
                Current_state = j-1;                                               % ex) current_path   = [0 0]
                Prev_state_zero = mod(2*Current_state, n_states) + 1;        % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(2*Current_state+1, n_states) + 1;      % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                % 입력 0을 받아서 현재 값이 되었을 경우
                if Current_state <  n_states / n_inputSymbols             
                    BM_zero = sum((output_zero(Prev_state_zero, :) - received_vector(t-1,:)).^2);  % 끝자리가 0인 코드에 0이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = sum((output_zero(Prev_state_one, :) - received_vector(t-1, :)).^2);  % 끝자리가 1인 코드에 0이 들어왔을 때
                    [Path_metric(j, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                             Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(j, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end                             

                % 입력 1을 받아서 현재 값이 되었을 경우    
                else                                                          
                    BM_zero = sum((output_one(Prev_state_zero, :) - received_vector(t-1,:)).^2);  % 끝자리가 0인 코드에 1이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = sum((output_one(Prev_state_one, :) - received_vector(t-1, :)).^2);  % 끝자리가 1인 코드에 1이 들어왔을 때
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
            % 그 이외에 나머지 stage
            for j = 1 : n_states
                % ex) current_path   = [0 0]
                Current_state = j-1;                                 
                Prev_state_zero = mod(2*Current_state, n_states)+1;          % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(2*Current_state+1, n_states)+1;        % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                % 입력 0을 받아서 현재 값이 되었을 경우
                if Current_state < n_states / n_inputSymbols                         
                    BM_zero = sum((output_zero(Prev_state_zero, :) - received_vector(t-1,:)).^2);  % 끝자리가 0인 코드에 0이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = sum((output_zero(Prev_state_one, :) - received_vector(t-1, :)).^2);  % 끝자리가 1인 코드에 0이 들어왔을 때
                    [Path_metric(j, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                       Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(j, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end

                % 입력 1을 받아서 현재 값이 되었을 경우
                else                            
                    BM_zero = sum((output_one(Prev_state_zero, :) - received_vector(t-1,:)).^2);  % 끝자리가 0인 코드에 1이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = sum((output_one(Prev_state_one, :) - received_vector(t-1, :)).^2);  % 끝자리가 1인 코드에 1이 들어왔을 때
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
    Survivor_path = Survivor_path(:, n_inputSymbols : end);
    Message_bit = Message_bit(:, n_inputSymbols : end);
    for t = num_message_bit + tail_bit :-1 : 1
        decoded_output(1,t) = Message_bit(current_state, t);
        current_state = Survivor_path(current_state, t);
    end
    decoded_output = decoded_output(1:num_message_bit);
end
