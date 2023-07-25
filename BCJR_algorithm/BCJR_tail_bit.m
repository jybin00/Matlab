%********************************************************************************%
%                                     BCJR_conv Decoder
% 
% This algorithm is reserved to the implementation of the Bahl, Cocke, Jelinek and Raviv (BCJR)
% algorithm. This function takes as input the channel output (corrupted data)
% and the a priori prob (we will set it to 1/2) and returns as output
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

    zero = repelem(-1, n);

    % Pri_p = 0.5; % The a priori probability.

    n_state = trellis.numStates;
    n_mem = log2(trellis.numStates);
    outputs = trellis.outputs;
    nextStates = trellis.nextStates;
    N0 = 2*sigma^2;

    @dec2bin_str;
    
    bin_outputs = 2 *int2bit(0:trellis.numOutputSymbols-1, n)' -1;
    
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
                x = bin_outputs(x+1, :);
                gamma(stage, prev_s, curr_s) = exp( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (N0) );
                input = input + 1;
            end 
            index = index + 1;
        end
        previous_s = reshape((nextstates + 1), 1, []);
    end

    for stage = n_mem + 1 : (N*R - n_mem)
        for prev_s = 1 : n_state
            current_s = nextStates(prev_s,:);
            input = 1;
            for curr_s = current_s + 1
                x = outputs(prev_s, input);
                x = bin_outputs(x+1, :);
                % using product sum instead of square.
                gamma(stage, prev_s, curr_s) = exp( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (N0) );
                input = input + 1;
            end
        end
    end

    previous_s = (1:2^n_mem);
    for stage = (N*R - n_mem) + 1 : N*R
        for prev_s = previous_s
            for curr_s = nextStates(prev_s, 0 + 1) + 1
                x = outputs(prev_s, 0 + 1);
                x = bin_outputs(x+1, :);
                % using product sum instead of square.
                gamma(stage, prev_s, curr_s) = exp( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (N0) );
            end 
        end
        previous_s = unique(nextStates(previous_s, 0 + 1))' + 1;
    end
    
    %************************************** alpha recursions********************************************
    alpha = zeros(N*R+1, n_state);
    alpha(1,1) = 1;

    for stage = 2 : N*R+1 % - n_mem
        for prev_s = 1 : n_state
            for curr_s = 1 : n_state
                alpha(stage, curr_s) = alpha(stage, curr_s) + gamma(stage-1, prev_s, curr_s) * alpha(stage-1, prev_s) ;
            end
        end
        alpha(stage,:) =  alpha(stage,:) / sum(alpha(stage,:)); % Normalization
    end

    % previous_s = (1:2^n_mem);
    % for stage = N*R+2 - n_mem : N*R + 1
    %     for prev_s = previous_s
    %         for curr_s = nextStates(prev_s, 0 + 1) + 1
    %             alpha(stage, curr_s) = alpha(stage, curr_s) + gamma(stage-1, prev_s, curr_s) * alpha(stage-1, prev_s) ;
    %         end 
    %     end
    %     previous_s = unique(nextStates(previous_s, 0 + 1))' + 1;
    %     alpha(stage,:) =  alpha(stage,:) / sum(alpha(stage,:)); % Normalization
    % end
    
    %************************************* beta recursions**********************************************
    beta = zeros(N*R+1, n_state);
    beta(N*R+1,1) = 1;
    j=0;
    % for stage = N*R+1 : -1 : N*R+1 - (n_mem)
    %     for prev_s = 1 : n_state/(n_mem-j)
    %         for curr_s = current_s
    %             beta(stage-1, prev_s) = beta(stage-1, prev_s) + gamma(stage-1, prev_s, curr_s) * beta(stage, curr_s);
    %         end
    %     end
    %     j = j + 1;
    % 
    %     beta(stage-1, :) = beta(stage-1, :) / sum(beta(stage-1,:)); 
    %     if(j >= n_mem) 
    %         break
    %     else 
    %         current_s = (1 : n_state/(n_mem-j));
    %     end
    % end
    % 
    % for stage = N*R+1 - (n_mem) :-1 : 2 + n_mem
    %     for prev_s = 1 : n_state
    %         for curr_s = 1 : n_state
    %             beta(stage-1, prev_s) = beta(stage-1, prev_s) + gamma(stage-1, prev_s, curr_s) * beta(stage, curr_s) ;
    %         end
    %     end
    %     beta(stage-1, :) = beta(stage-1, :) / sum(beta(stage-1,:)) ; % Normalization
    % end
    % 
    % for stage = 1 + n_mem : -1 : 2
    %     for prev_s = 1: n_state
    %         for curr_s = 1 : n_state
    %             beta(stage-1, prev_s) = beta(stage-1, prev_s) + gamma(stage-1, prev_s, curr_s) * beta(stage, curr_s);
    %         end
    %     end
    %     beta(stage-1, :) = beta(stage-1, :) / sum(beta(stage-1,:)) ; % Normalization
    % end
    
    for stage = N/n +1 :-1 : 2
        for prev_s = 1 : n_state
            for curr_s = 1 : n_state
                beta(stage-1, prev_s) = beta(stage-1, prev_s) + gamma(stage-1, prev_s, curr_s) * beta(stage, curr_s) ;
            end
        end
        beta(stage-1, :) = beta(stage-1, :) / sum(beta(stage-1,:)) ; % Normalization
    end
    %------------------------------------- Computing the LLRs -------------------------------------%
    for stage= 1:N*R
        up= 0;
        down= 0;
        for prev_s = 1 : n_state
            current_s = nextStates(prev_s,:) + 1;
            in = 1;
            for curr_s = current_s
                if     (in == 2) % input = 1
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
            x = zero;
        else
            tmp_x = dec2bin(codeword, n);
            bin = zeros(1, n);
            m = 1;
            while m < n + 1
                bin(m) = 2*(tmp_x(m) > '0') - 1; 
                m = m + 1; 
            end
            x = bin;
        end
    end
end
