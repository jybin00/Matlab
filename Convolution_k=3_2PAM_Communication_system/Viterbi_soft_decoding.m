% Viterbi decoding
% Convolution code는 다음과 같음. (3, [7, 5]) = (7, [111], [101])

function decoded_output = Viterbi_soft_decoding(received_signal, num_message_bit, Es)  
    
    tail_bit = 2;
    decoded_output = zeros(1, num_message_bit + tail_bit);  % demodulated_bit + 2bit(tail bits)
%---------------------------------------------------------------------------------
    % branch metric 계산을 위해서 각 state가 0과 1이 입력으로 들어왔을 때 어떤 출력을 하는지 미리 계산
    output_zero =  [-sqrt(Es) -sqrt(Es);  sqrt(Es) sqrt(Es);
                    sqrt(Es) -sqrt(Es);  -sqrt(Es) sqrt(Es)];
     
    % 1이 입력으로 들어왔을 때 p1, p2

    output_one =  [sqrt(Es) sqrt(Es);  -sqrt(Es) -sqrt(Es);	
                   -sqrt(Es) sqrt(Es);	sqrt(Es) -sqrt(Es)]; 
    
    codeword = received_signal;
%--------------------------------------------------------------------------------------
    Path_metric = inf(4, num_message_bit + tail_bit + 1);                         % 4 X num_message_bit + 2 + 1 Path metric 저장 위한 matrix
    Message_bit = inf(4, num_message_bit + tail_bit + 1);   % path message array
    Survivor_path = zeros(4, num_message_bit + tail_bit + 1); 
    for t = 1 : num_message_bit+ tail_bit + 1                                            % tail bits 존재
        if t == 1
            Path_metric(1, 1) = 0;	                                                     % 맨 처음 state
        elseif t > 1 && t <= tail_bit                                                           % t = 2는 state가 펼쳐지는 시간 
            for j = 1 : 2^(t-1)
                Current_state = (4/2^(t-1))*(j-1);                                 % ex) current_path   = [0 0], [1 0]
                Prev_state_zero = mod(2*Current_state, 4) + 1;              % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(2*Current_state+1, 4) +1;       % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_state < 2    % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = sum((output_zero(Prev_state_zero, :) - codeword(t-1,:)).^2);  % 끝자리가 0인 코드에 0이 들어왔을 때 outputd으로 BM 계산
                    BM_one  = sum((output_zero(Prev_state_one, :) - codeword(t-1, :)).^2);  % 끝자리가 1인 코드에 0이 들어왔을 때
                    [Path_metric(1+Current_state, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                                Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(1+Current_state, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end   
				else	                                                                             % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = sum((output_one(Prev_state_zero, :) - codeword(t-1,:)).^2);
                    BM_one  = sum((output_one(Prev_state_one, :) - codeword(t-1, :)).^2);
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
            for j = 1 : 2^(num_message_bit + 2 + 1 -t)			
                Current_state = j-1;                                                      % ex) current_path   = [0 0]
                Prev_state_zero = mod(2*Current_state, 4)+1;             % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(2*Current_state + 1, 4) +1;      % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_state < 2                                                    % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = sum((output_zero(Prev_state_zero, :) - codeword(t-1,:)).^2);
                    BM_one  = sum((output_zero(Prev_state_one, :) - codeword(t-1, :)).^2);
                    [Path_metric(j, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                             Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(j, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end         
                else                                                                             % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = sum((output_one(Prev_state_zero, :) - codeword(t-1,:)).^2);
                    BM_one  = sum((output_one(Prev_state_one, :) - codeword(t-1, :)).^2);
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
            for j = 1 : 4
                Current_state = j-1;                                                        % ex) current_path   = [0 0]
                Prev_state_zero = mod(2*Current_state, 4)+1;               % prev_path_from 0 = [0 0] 끝자리 0
                Prev_state_one  = mod(2*Current_state + 1, 4) +1;        % prev_path_from 1 = [0 1] 끝자리 1

                % BM = min{ 끝자리 0에서 온 PM + BM, 끝자리 1에서 온 PM + BM}
                if Current_state < 2    % 입력 0을 받아서 현재 값이 되었는가?
                    BM_zero = sum((output_zero(Prev_state_zero, :) - codeword(t-1,:)).^2);
                    BM_one  = sum((output_zero(Prev_state_one, :) - codeword(t-1, :)).^2);
                    [Path_metric(j, t), I] = min([Path_metric(Prev_state_zero, t-1) + BM_zero, ...
                                                       Path_metric(Prev_state_one, t-1) + BM_one]);
                    Message_bit(j, t) = 0;
                    if I > 1
                        Survivor_path(1+Current_state, t) = Prev_state_one;                        
                    else
                        Survivor_path(1+Current_state, t) = Prev_state_zero;
                    end
                else                            % 입력 1을 받아서 현재 값이 되었는가?
                    BM_zero = sum((output_one(Prev_state_zero, :) - codeword(t-1,:)).^2);
                    BM_one  = sum((output_one(Prev_state_one, :) - codeword(t-1, :)).^2);
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
    current_state = 00 + 1;                         % 맨 마지막 state = [0 0]
    Survivor_path = Survivor_path(:,2:end);
    Message_bit = Message_bit(:, 2:end);
    for t = num_message_bit + tail_bit :-1 : 1
        decoded_output(1,t) = Message_bit(current_state, t);
        current_state = Survivor_path(current_state, t);
    end
    decoded_output = decoded_output(1:num_message_bit);
end
