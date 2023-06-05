function encoder = layered_encoder_w_transpose(N, K, eps, input_type)
%% Upper triangular transform이 되도록.
L = log2(N);

layer_A = cell(1, L);
layer_C = cell(1, L);
msg_len = zeros(1, L);

%% Parameter setting
switch N
    case 16
        switch K
            case 8
                %% (16,8)
                polarization = [false, false, true, true];
%                 num_info_list = [0, 0, 5, 3];
%                 num_conn_list = [0, 0, 0, 8];     
                layer_A{3} = [1,2,3,4,5];
                layer_A{4} = [14,15,16];
                layer_C{4} = [4,6,7,8,10,11,12,13];
            otherwise
                error('Not implemented encoder configuration: (N, K) = (%d, %d)', N, K);
        end
    case 32
        switch K
            case 8
                %% (32,8)
                polarization = [false, false, true, false, true];
%                 num_info_list = [0, 0, 2, 0, 6];
%                 num_conn_list = [0, 0, 0, 0, 8];
                layer_A{3} = [1,2];
                layer_A{5} = [16,24,28,30,31,32];
                layer_C{5} = [14,15,20,22,23,26,27,29];
            case 12
                %% (32,12)
                polarization = [false, false, true, false, true];
%                 num_info_list = [0, 0, 4, 0, 8];
%                 num_conn_list = [0, 0, 0, 0, 8];
                layer_A{3} = [1,2,3,5];
                layer_A{5} = [16,24,27,28,29,30,31,32];
                layer_C{5} = [8,12,14,15,20,22,23,26];
            case 16
                %% (32,16)
                polarization = [false, true, true, true, true];
%                 num_info_list = [0, 2, 2, 4, 8];
%                 num_conn_list = [0, 0, 4, 8, 16];
            otherwise
                error('Not implemented encoder configuration: (N, K) = (%d, %d)', N, K);
        end
    case 64
        switch K
            case 8
                %% (64,8)
                polarization = [false, false, false, true, false, true];
%                 num_info_list = [0, 0, 0, 2, 0, 6];
%                 num_conn_list = [0, 0, 0, 0, 0, 16];
                layer_A{4} = [1,2];
                layer_A{6} = [48,56,60,62,63,64];
                layer_C{6} = [32,  16,24,28,30,31,40,44,46,47,52,54,55,58,59,61];
            case 12
                %% (64,12)
                polarization = [false, true, false, true, false, true];
%                 num_info_list = [0, 1, 0, 5, 0, 6];
%                 num_conn_list = [0, 0, 0, 4, 0, 16];
                layer_A{2} = [1];
                layer_A{4} = [1,2,3,5,9];
                layer_A{6} = [48,56,60,62,63,64];
                layer_C{4} = [7,10,11,13];
                layer_C{6} = [32,  16,24,28,30,31,40,44,46,47,52,54,55,58,59,61];
            case 16
                %% (64,16)
                polarization = [false, false, true, false, false, true];
%                 num_info_list = [0, 0, 4, 0, 0, 12];
%                 num_conn_list = [0, 0, 0, 0, 0, 8];
                layer_A{3} = [1,2,3,5];
                layer_A{6} = [32,48,54,55,56,58,59,60,61,62,63,64];
                layer_C{6} = [28,30,31,40,44,46,47,52];
            otherwise
                error('Not implemented encoder configuration: (N, K) = (%d, %d)', N, K);
        end
    case 128
        switch K
            case 12
                %% (128,12)
                polarization = [false, false, false, true, false, false, true];
%                 num_info_list = [0, 0, 0, 4, 0, 0, 8];
%                 num_conn_list = [0, 0, 0, 0, 0, 0, 16];
                layer_A{4} = [1,2,3,5];
                layer_A{7} = [64,96,112,120,124,126,127,128];
                layer_C{7} = [62,63,88,92,94,95,104,108,110,111,116,118,119,122,123,125];
            case 16
                %% (128,16)
                polarization = [false, false, false, true, false, false, true];
%                 num_info_list = [0, 0, 0, 5, 0, 0, 11];
%                 num_conn_list = [0, 0, 0, 0, 0, 0, 16];
                layer_A{4} = [1,2,3,5,9];
                layer_A{7} = [64 96 112 120 124 126 127 128, 122, 123, 125];
                layer_C{7} = [56 60 62 63 80 88 92 94 95 104 108 110 111 116 118 119];  
            otherwise
                error('Not implemented encoder configuration: (N, K) = (%d, %d)', N, K);
        end
    otherwise
        error('Not implemented encoder configuration: (N, K) = (%d, %d)', N, K);
end

%% Encoder setting
G = cell(1, L);
% layer_A = cell(1, L);
% layer_C = cell(1, L);
% msg_len = nan(1, L);
for idx_l_outer = 1:L
%     if idx_l_outer == L
%         is_transpose = false;
%     else
%         is_transpose = true;
%     end
    G{idx_l_outer} = polar.generator_matrix(2^idx_l_outer);
%     [~, layer_C{idx_l_outer}, layer_A{idx_l_outer}] = ...
%         polar_partition(eps, 2^idx_l_outer, num_conn_list(idx_l_outer), num_info_list(idx_l_outer), is_transpose);
    msg_len(idx_l_outer) = length(layer_A{idx_l_outer});
end

if strcmpi(input_type, 'int')
    encoder = @(msg_in) encode_layer(int2bit(msg_in, K).');
elseif strcmpi(input_type, 'bit')
    encoder = @encode_layer;
else
    error("Possible input type = 'int' or 'bit', current: %s.", input_type)
end

%% Master Encoder
    function codeword = encode_layer(msg_in)
        msg_l_tmp = cell(1, L);
        for idx_l = 1:L
            msg_l_tmp{idx_l} = msg_in(1:msg_len(idx_l));
            msg_in = msg_in(msg_len(idx_l)+1:end);
        end

        assert(isempty(msg_in));

        msg_l = cell(1, L);
        msg_l_enc = cell(1, L);
        connecting_msg_tmp = [];
        for idx_l = 1:L
            if polarization(idx_l)
                msg_l{idx_l} = zeros(1, 2^idx_l);
                msg_l{idx_l}(layer_C{idx_l}) = connecting_msg_tmp;
                msg_l{idx_l}(layer_A{idx_l}) = msg_l_tmp{idx_l};
                if idx_l == L
                    msg_l_enc{idx_l} = mod(msg_l{idx_l} * G{idx_l}, 2);
                else
                    msg_l_enc{idx_l} = mod(msg_l{idx_l} * G{idx_l}.', 2);
                end
            else
                msg_l_enc{idx_l} = connecting_msg_tmp;
            end
            connecting_msg_tmp = msg_l_enc{idx_l};
        end

        codeword = msg_l_enc{L};
    end

%% ============ Utility function ============
    function [frozen_idx, con_idx, info_idx] = polar_partition(eps, n, con_num, info_num, is_transpose)
        cap_sub = 1 - arrayfun(@(x) polar.BEC_transformed_equivalent(eps, x-1, log2(n)), 1:n);
        [~, idx] = sort(cap_sub, 'ascend');
        idx = flip(idx);

        if is_transpose
            frozen_num = n - info_num - con_num;
            frozen_idx = idx(1:frozen_num);
            con_idx = idx(frozen_num+1:frozen_num+con_num);
            info_idx = idx(frozen_num+con_num+1:end);
        else
            info_idx = idx(1:info_num);
            con_idx = idx(info_num+1:info_num+con_num);
            frozen_idx = idx(info_num+con_num+1:end);
        end

        info_idx = sort(info_idx);
        con_idx = sort(con_idx);
        frozen_idx = sort(frozen_idx);
    end
end