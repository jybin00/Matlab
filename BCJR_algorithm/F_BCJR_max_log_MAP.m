

function LLR = F_BCJR_max_log_MAP(y, trellis, sigma)

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
    N0 = 2*sigma^2;

    bin_outputs = 2 *int2bit(0:trellis.numOutputSymbols-1, n)' -1;
    
    %********* Computing gamma for all states at each time ***************
    Gamma = -inf(N*R, n_state, n_state); % we suppose that the first state is the 0 state which can be handled at the encoders.
    previous_s = 1;
    for stage = 1 : n_mem
        nextstates = nextStates(previous_s, :);
        index = 1;
        for prev_s = previous_s
            input = 1;
            for curr_s = nextstates(index, :) + 1
                x = outputs(prev_s, input);
                x = bin_outputs(x+1, :);
                Gamma(stage, prev_s, curr_s) = ( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (N0) );
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
                Gamma(stage, prev_s, curr_s) = ( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (N0) );
                input = input + 1;
            end
        end
    end

    previous_s = (1:2^n_mem);
    index_k = 0;
    for stage = (N*R - n_mem) + 1 : N*R
        for prev_s = previous_s
            for curr_s = nextStates(prev_s, 0 + 1) + 1
                x = outputs(prev_s, 0 + 1);
                x = bin_outputs(x+1, :);
                % using product sum instead of square.
                Gamma(stage, prev_s, curr_s) = ( 4 * sum(x .* y(n*(stage-1)+1 : n*stage)) / (N0) );
            end 
        end
        index_k = index_k + 1;
        previous_s = (1:2^(n_mem-index_k));
    end
    
    %************** alpha recursions********************

    Alpha = -inf(N*R+1, n_state);
    Alpha(1,1) = 0;
    for stage = 2 : N*R+1
        for curr_s = 1 : n_state
            for prev_s = 1 : n_state
                Alpha(stage, curr_s) = max(Alpha(stage, curr_s), (Gamma(stage-1, prev_s, curr_s) + Alpha(stage-1, prev_s))) ;
            end
        end
        %Alpha(stage,:) =  max(Alpha(stage,:)); 
    end
    
    %************** beta recursions********************
    Beta = -inf(N*R+1, n_state);
    Beta(N*R+1,:) = 0;

    for stage = N/n +1 :-1 : 2
        for prev_s = 1 : n_state
            current_s = nextStates(prev_s,:) + 1;
            for curr_s = current_s
                Beta(stage-1, prev_s) = max(Beta(stage-1, prev_s), (Gamma(stage-1, prev_s, curr_s) + Beta(stage, curr_s))) ;
            end
        end
        %Beta(stage-1, :) = max(Beta(stage-1, :)); 
    end
    
    %%%%%%%%%%%%%%%%%%% Computing the LLRs %%%%%%%%%%%%%%%
    for stage= 1:N*R
        up= 0;
        down= 0;
        for prev_s = 1 : n_state
            current_s = nextStates(prev_s,:) + 1;
            in = 1;
            for curr_s = current_s
                if     (in == 2) % input = 1
                   up = max(up, Alpha(stage, prev_s) + Gamma(stage, prev_s, curr_s) + Beta(stage+1, curr_s));
                elseif (in == 1) % input = 0
                   down = max(down, Alpha(stage, prev_s) + Gamma(stage, prev_s, curr_s) + Beta(stage+1, curr_s));
                end
                in = in + 1;
            end
        end
        LLR(stage) = up - down;
    end
end
