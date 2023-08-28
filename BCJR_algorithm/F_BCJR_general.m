%% BCJR general decoder
% trellis 정보 받아서 preprocessing하고 decoder를 반환함.


function decoder = F_BCJR_general(trellis, sigma, varargin)

    n_state = trellis.numStates;
    n_mem = log2(trellis.numStates);
    outputs = trellis.outputs;
    nextStates = trellis.nextStates;
    N0 = 2*sigma^2;

    n = log2(trellis.numOutputSymbols); % n = 2
    k = log2(trellis.numInputSymbols); % k=1
    R = k/n; % coding rate, R=1/2

    bin_outputs = 2 *int2bit(0:trellis.numOutputSymbols-1, n)' -1;

    %%

    p = inputParser;

    default_type = 'True BCJR';
    valid_type = {'True BCJR', 'log-max-MAP', 'log-MAP', 'APP BCJR', 'APP max-log', 'APP log'};
    check_type = @(x) any(validatestring(x, valid_type));

    addRequired(p, 'trellis', @istrellis);
    addRequired(p, 'sigma', @isnumeric);
    addOptional(p, 'type', default_type, check_type);

    parse(p, trellis, sigma, varargin{:});

    if strcmp(p.Results.type, 'True BCJR')
        decoder = @BCJR_tail_bit;

    elseif strcmp(p.Results.type, 'log-max-MAP')
        decoder = @BCJR_max_log_MAP;

    elseif strcmp(p.Results.type, 'log-MAP')
        decoder = @BCJR_log_MAP;

    elseif strcmp(p.Results.type, 'APP BCJR')
        decoder = comm.APPDecoder(...
            'TrellisStructure', trellis, ...
            'Algorithm', 'Max', ...
            'CodedBitLLROutputPort', false, ...
            'TerminationMethod','Terminated');

    elseif strcmp(p.Results.type, 'APP max-log')
        decoder = comm.APPDecoder(...
            'TrellisStructure', trellis, ...
            'Algorithm', 'Max', ...
            'CodedBitLLROutputPort', false, ...
            'TerminationMethod','Terminated');
        
    elseif strcmp(p.Results.type, 'APP log')
        decoder = comm.APPDecoder(...
            'TrellisStructure', trellis, ...
            'Algorithm', 'Max*', ...
            'CodedBitLLROutputPort', false, ...
            'TerminationMethod','Terminated');
    end

    %%
    function LLR = BCJR_tail_bit(size, y)
        
        tmp_size = size;
        y = y';

        N = length(y); % y is the channel output

        LLR = zeros(1, N*R);
    
        zero = repelem(-1, n);
    
        % Pri_p = 0.5; % The a priori probability.
        
        
        
        %*************************** Computing gamma for all states at each time *****************************
        gamma = zeros(N*R, n_state, n_state); % we suppose that the first state is the 0 state which can be handled at the encoders.
        coder.varsize('previous_s');
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
            previous_s = reshape((nextstates + 1), 1, 2^stage);
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
        index_k = 0;
        for stage = (N*R - n_mem) + 1 : N*R
            for prev_s = previous_s
                for curr_s = nextStates(prev_s, 0 + 1) + 1
                    x = outputs(prev_s, 0 + 1);
                    x = bin_outputs(x+1, :);
                    % using product sum instead of square.
                    gamma(stage, prev_s, curr_s) = exp( 4*sum(x .* y(n*(stage-1)+1 : n*stage) ) / (N0) );
                end 
            end
            index_k = index_k + 1;
            previous_s = (1:2^(n_mem-index_k));
        end
        
        %************************************** alpha recursions********************************************
        alpha = zeros(N*R+1, n_state);
        alpha(1,1) = 1;
    
        for stage = 2 : N*R+1 % - n_mem
            for prev_s = 1 : n_state
                current_s = nextStates(prev_s,:);
                for curr_s = current_s + 1
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
                current_s = nextStates(prev_s,:);
                for curr_s = current_s + 1
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
    
    %%
    function LLR = BCJR_max_log_MAP(size, y)

        tmp_size = size;
        y = y';
    
        N = length(y); % y is the channel output
        LLR = zeros(1, N*R);
    
        % Pri_p = 0.5; % The a priori probability.
        
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

end

