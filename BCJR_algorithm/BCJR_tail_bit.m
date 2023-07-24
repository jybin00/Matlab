%********************************************************************************%
%                                     BCJR_conv Decoder
% 
% This algorithm is reserved to the implementation of the Bahl, Cocke, Jelinek and Raviv (BCJR)
% algorithm. This function takes as input the channel output (corrupted
% data) and the a priori prob (we will set it to 1/2) and returns as output
% the APP Log Liklihood Ratio (LLR) for every data input. It is usually called a
% Soft Input Soft Output (SISO) decoder. It can be applied to any code
% having a finite state machine, in our case we will use it for rate-1/n convolutional codes.
%         					
%                                              K. Elkhalil, SUP'COM Tunisia
%                				             
% *******************************************************************************%

function LLR = BCJR_tail_bit(y, trellis, sigma)

    N = length(y); % y is the channel output
    n = log2(trellis.numOutputSymbols); % n = 2
    k = log2(trellis.numInputSymbols); % k=1
    R = k/n; % coding rate, R=1/2
    LLR = zeros(1, N*R);

    % Pri_p = 0.5; % The a priori probability.

    n_state = trellis.numStates;
    n_mem = log2(trellis.numStates);
    outputs = trellis.outputs;
    nextStates = trellis.nextStates;
    
    opening = cell(1,n_mem);
    index = 1;
    for i = 1:n_mem
        for j = index+1
            opening{i} = [opening{i}; trellis.nextStates(j,:)];
        end
        index = cell2mat(opening(i));
    end

    @dec2bin_str;
    
    %*************************** Computing gamma for all states at each time *****************************
    gamma = zeros(N*R, n_state, n_state); % we suppose that the first state is the 0 state which can be handled at the encoders.
    previous_s = 1;
    for stage = 1 : n_mem
        nextstates = nextStates(previous_s, :);
        index = 1;
        for prev_s = previous_s
            input = 1;
            for curr_s = nextstates(index, :) + 1
                x = outputs(prev_s, input);
                x = dec2bin_str(x);
                gamma(stage, prev_s, curr_s) = exp( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (2*sigma^2) );
                input = input + 1;
            end 
            index = index + 1;
        end
        previous_s = nextstates + 1;
    end
    for stage = n_mem + 1 : N*R - n_mem
        for prev_s = 1 : n_state
            current_s = nextStates(prev_s,:);
            input = 1;
            for curr_s = current_s + 1
                x = outputs(prev_s, input);
                x = dec2bin_str(x);
                % using product sum instead of square.
                gamma(stage, prev_s, curr_s) = exp( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (2*sigma^2) );
                input = input + 1;
            end
        end
    end

    previous_s = (1:2^n_mem);
    for stage = N*R - n_mem + 1 : N*R
        for prev_s = previous_s
            for curr_s = nextStates(prev_s, 0 + 1) + 1
                x = outputs(prev_s, 0 + 1);
                x = dec2bin_str(x);
                % using product sum instead of square.
                gamma(stage, prev_s, curr_s) = exp( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (2*sigma^2) );
            end 
        end
        previous_s = unique(nextStates(previous_s, 0 + 1))' + 1;
    end
    
    %************************************** alpha recursions********************************************
    alpha = zeros(N*R+1, n_state);
    alpha(1,1) = 1;

    for stage = 2 : N*R+1 - n_mem
        for prev_s = 1 : n_state
            for curr_s = 1 : n_state
                alpha(stage, curr_s) = alpha(stage, curr_s) + gamma(stage-1, prev_s, curr_s) * alpha(stage-1, prev_s) ;
            end
        end
        alpha(stage,:) =  alpha(stage,:) / sum(alpha(stage,:)); % Normalization
    end

    previous_s = (1:2^n_mem);
    for stage = N*R+2 - n_mem : N*R + 1
        for prev_s = previous_s
            for curr_s = nextStates(prev_s, 0 + 1) + 1
                alpha(stage, curr_s) = alpha(stage, curr_s) + gamma(stage-1, prev_s, curr_s) * alpha(stage-1, prev_s) ;
            end 
        end
        previous_s = unique(nextStates(previous_s, 0 + 1))' + 1;
        alpha(stage,:) =  alpha(stage,:) / sum(alpha(stage,:)); % Normalization
    end
    
    %************************************* beta recursions**********************************************
    beta = zeros(N*R+1, n_state);
    beta(N*R+1,1) = 1;
    j=0;  current_s = 1;
    for stage = N*R+1 : -1 : N*R+1 - (n_mem)
        for prev_s = 1 : n_state/(n_mem-j)
            for curr_s = current_s
                beta(stage-1, prev_s) = beta(stage-1, prev_s) + gamma(stage-1, prev_s, curr_s) * beta(stage, curr_s);
            end
        end
        j = j + 1;

        if(j >= n_mem) 
            break
        else 
            current_s = (1 : n_state/(n_mem-j));
        end
        beta(stage-1, :) = beta(stage-1, :) / sum(beta(stage-1,:)); 
    end

    for stage = N*R+1 - (n_mem) :-1 : 2 + n_mem
        for prev_s = 1 : n_state
            for curr_s = 1 : n_state
                beta(stage-1, prev_s) = beta(stage-1, prev_s) + gamma(stage-1, prev_s, curr_s) * beta(stage, curr_s) ;
            end
        end
        beta(stage-1, :) = beta(stage-1, :) / sum(beta(stage-1,:)) ; % Normalization
    end

    for stage = 1 + n_mem : -1 : 2
        for prev_s = 1: n_state
            for curr_s = 1 : n_state
                beta(stage-1, prev_s) = beta(stage-1, prev_s) + gamma(stage-1, prev_s, curr_s) * beta(stage, curr_s);
            end
        end
        beta(stage-1, :) = beta(stage-1, :) / sum(beta(stage-1,:)) ; % Normalization
    end
    
    %%%%%%%%%%%%%%%%%%% Computing the LLRs %%%%%%%%%%%%%%%
    for stage= 1:N*R
        up= 0;
        down= 0;
        for prev_s = 1 : n_state
            current_s = nextStates(prev_s,:);
            in = 1;
            for i = current_s
                curr_s = i + 1;
                if (in == 2) % input = 1
                   up = up + alpha(stage, prev_s) * gamma(stage, prev_s, curr_s) * beta(stage+1, curr_s);
                elseif (in == 1) % input = 0
                   down = down + alpha(stage, prev_s) * gamma(stage, prev_s, curr_s) * beta(stage+1, curr_s);
                end
                in = in + 1;
            end
        end
        LLR(stage) = log(up/down);
    end

    function x = dec2bin_str(codeword)
        if codeword == 0
            x = repelem(-1, n);
        else
            tmp_x = dec2bin(codeword, n);
            bin = zeros(1, n);
            i = 1;
            while i < n + 1
                bin(i) = 2*(tmp_x(i) > '0') - 1; 
                i = i + 1; 
            end
            x = bin;
        end
    end
end