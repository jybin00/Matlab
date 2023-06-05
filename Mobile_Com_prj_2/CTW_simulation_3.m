clc; clear; close all;

TRIAL = 1e7;
GATHERING = 1e3;

%% Parameters of polar codes
N = 128;
K = 16;
channel_mode = 'BEC';
switch channel_mode
    case 'BEC'
        eps_list = 0.85:-0.1:0.65;
        eps_list = 0.7;
    otherwise
        error('mode wrong')
end

msg_enum_ML = (1:2^K)-1;
msg_enum_max = msg_enum_ML(end);

%% For each BSC channel
FER_polar_ML = nan(size(eps_list));
FER_layer_ML = nan(size(eps_list));
for idx_eps = 1:length(eps_list)
    %% Frozen bits selection
    eps = eps_list(idx_eps);

    switch channel_mode
        case 'BEC'
            [frozen_idx, info_idx] = polar.frozen_bit_selection(N, K, eps, 'BEC');
%             info_idx = [64 96 112 120 124 126 127 128, 110, 111, 116, 118, 119, 122, 123, 125];
            eps_enc = eps;
        otherwise
            error('mode wrong')
    end

    G = polar.generator_matrix(N);
    G = G(info_idx, :);

    %% Dictionary generation
    encode_polar = @(msg_enum) mod(int2bit(msg_enum, K).' * G, 2);
    codeword_dict_polar = encode_polar(msg_enum_ML);
    decode_polar_ML = ML_decoder(codeword_dict_polar, channel_mode);

    encode_layer = layered_encoder_w_transpose(N,K,eps_enc,'int');
    codeword_dict_layer = nan(size(codeword_dict_polar));
    for idx_dict = 1:size(codeword_dict_layer, 1)
        codeword_dict_layer(idx_dict, :) = encode_layer(idx_dict-1);
    end
    decode_layer_ML = ML_decoder(codeword_dict_layer, channel_mode);

    %% Main loop
    switch channel_mode
        case 'BEC'
            channel = @(x) binary_erasure_channel(eps, x);
        otherwise
            error('mode wrong')
    end

    trial = 0;
%     error_num_polar_ML = 0;
    error_num_layer_ML = 0;
    while trial < TRIAL && error_num_layer_ML < GATHERING
        tic;
        foreachtrial = 1e2;
        msg_enum = randi([0, msg_enum_max], 1, foreachtrial);
%         dec_msg_enum_polar_ML = arrayfun(@(msg_enum) decode_polar_ML(channel(encode_polar(msg_enum))), msg_enum);
        dec_msg_enum_layer_ML = arrayfun(@(msg_enum) decode_layer_ML(channel(encode_layer(msg_enum))), msg_enum);

        trial = trial + foreachtrial;
%         error_num_polar_ML = error_num_polar_ML + sum(~(dec_msg_enum_polar_ML == msg_enum));
        error_num_layer_ML = error_num_layer_ML + sum(~(dec_msg_enum_layer_ML == msg_enum));

        toc;
%         fprintf('===\n%s %.2f [%d / %d]\nFrame %d / %d \nError (polar ML, layer ML) (%d, %d) / %d\n', ...
%             channel_mode, eps, idx_eps, length(eps_list), trial, TRIAL, error_num_polar_ML, error_num_layer_ML, GATHERING);
        fprintf('===\n%s %.2f [%d / %d]\nFrame %d / %d \nError (layer ML) (%d) / %d\n', ...
            channel_mode, eps, idx_eps, length(eps_list), trial, TRIAL, error_num_layer_ML, GATHERING);        
    end

%     FER_polar_ML(idx_eps) = error_num_polar_ML / trial;
    FER_layer_ML(idx_eps) = error_num_layer_ML / trial;

    save_path = "(N,K)=(" + string(N) + "," + string(K) + ") FER.mat";
    save(save_path, 'eps_list', 'FER_polar_ML', 'FER_layer_ML')

%     if FER_polar_ML(idx_eps) == 0 && FER_layer_ML(idx_eps) == 0
%         break;
%     end
end
