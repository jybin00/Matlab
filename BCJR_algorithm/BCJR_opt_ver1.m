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

function LLR = BCJR_opt_ver1(y, trellis, sigma)

    N = length(y); % y is the channel output
    n = log2(trellis.numOutputSymbols); % n = 2
    k = log2(trellis.numInputSymbols); % k=1
    R = k/n; % coding rate, R=1/2
    LLR = zeros(1, N*R);
    Pri_p = 0.5; % The a priori probability.
    n_state = trellis.numStates;
    
    %********* Computing gamma for all states at each time ***************
    gamma = zeros(N*R, n_state, n_state); % we suppose that the first state is the 0 state which can be handled at the encoders.
    for stage = 1 : N*R
        for prev_s = 1 : n_state
            current_s = trellis.nextStates(prev_s,:);
            input = 1;
            for in = current_s
                x = trellis.outputs(prev_s, input);
                if x == 0
                    x = repelem(-1, n);
                else
                    tmp_x = dec2bin(x, n);
                    bin = zeros(1, n);
                    i = 1;
                    while i < n + 1
                        bin(i) = 2*(tmp_x(i) > '0') - 1; 
                        i = i + 1; 
                    end
                    x = bin;
                end
                % using product sum instead of square.
                curr_s = in + 1;
                gamma(stage, prev_s, curr_s) = exp( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (2*sigma^2) );
                input = input + 1;
            end
        end
    end
    
    %************** alpha recursions********************
    alpha = zeros(N*R+1, trellis.numStates);
    alpha(1,1) = 1;
    for stage = 2 : N*R+1
        for prev_s = 1 : n_state
            for curr_s = 1 : n_state
                alpha(stage, curr_s) = alpha(stage, curr_s) + gamma(stage-1, prev_s, curr_s) * alpha(stage-1, prev_s) ;
            end
        end
        alpha(stage,:) =  alpha(stage,:) / sum(alpha(stage,:)); % Normalization
    end
    
    %************** beta recursions********************
    beta = zeros(N*R+1, n_state);
    beta(N*R+1,:) = alpha(N*R+1,:);
    for stage = N*R+1:-1:2
        for prev_s = 1 : n_state
            for curr_s = 1 : n_state
                beta(stage-1, prev_s) = beta(stage-1, prev_s) + gamma(stage-1, prev_s, curr_s) * beta(stage, curr_s) ;
            end
        end
        beta(stage-1, :) = beta(stage-1, :) / sum(beta(stage-1,:)) ; % Normalization
    end
    
    %%%%%%%%%%%%%%%%%%% Computing the LLRs %%%%%%%%%%%%%%%
    for stage= 1:N*R
        up= 0;
        down= 0;
        for prev_s = 1 : n_state
            current_s = trellis.nextStates(prev_s,:);
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
end