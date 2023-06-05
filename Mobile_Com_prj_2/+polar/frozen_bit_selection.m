function [frozen_idx, info_idx] = frozen_bit_selection(N,K,eps,channel)
% ============================================================
% BEC에 대해 frozen bit 위치 반환
% Input:
%   N: Codeword 길이
%   K: Message 길이
%      i.e., code rate = K / N
%   eps: BSC epsilon, BEC epsilon, AWGN snr
%   channel: 'BEC', 'BSC', 'AWGN'
% Output:
%   frozen_idx: frozen bit index
%   info_idx  : information bit index
% ============================================================

assert (N>0 && N>=K)
assert (K>0)
% assert (eps>=0 && eps<=1)

switch channel
    case 'BEC'
        capacity_sub = 1 - ...
            arrayfun(@(x) polar.BEC_transformed_equivalent(eps, x-1, log2(N)), 1:N);
    case 'BSC'
        equiv_eps = 2*sqrt(eps*(1-eps));
        capacity_sub = 1 - ...
            arrayfun(@(x) polar.BEC_transformed_equivalent(equiv_eps, x-1, log2(N)), 1:N);
    case 'AWGN'
        equiv_eps = qfunc(eps);
        capacity_sub = 1 - ...
            arrayfun(@(x) polar.AWGN_transformed_equivalent(equiv_eps, x-1, log2(N)), 1:N);
    otherwise
        error("Parameter channel '%s' is unexpected. Available only 'BEC' or 'BSC'.", channel)
end

[~, idx] = sort(capacity_sub, 'ascend');

frozen_idx = idx(1:(N-K));
info_idx = idx(N-K+1:end);
frozen_idx = sort(frozen_idx);
info_idx = sort(info_idx);

if strcmpi(channel, 'AWGN')
    [info_idx, frozen_idx] = Find_InformationSet_5GNR(N,K);
end

end