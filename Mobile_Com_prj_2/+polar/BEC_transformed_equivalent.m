function y = BEC_transformed_equivalent(epsilon, idx, N)
% ============================================================
% Binary erasure channel에서 polar transform에 의해 생성되는
% subchannel은 equivalent BEC인데 그것의 erasure probability 반환
% Input:
%   epsilon: BEC erasure probability
%   idx: channel index
%   N: channel 개수 = 2^N
% Output:
%   y: subchannel equiv. BSC prob
% Note:
%   For BSC, epsilon = 2 * sqrt(eps * (1-eps)) <- Bhattachayya parameter
% ============================================================

idx_bin = de2bi(idx, N);

assert (~isempty(idx_bin))
assert (N > 0)


idx_current = idx_bin(1);
% idx_current = idx_bin(end);
if length(idx_bin) == 1
    if idx_current == 0
        % bad channel
        y = 2 * epsilon - epsilon^2;
    else
        % good channel
        y = epsilon^2;
    end
    return
end

idx_remain = bi2de(idx_bin(2:end));
% idx_remain = bi2de(idx_bin(1:end-1));
y = polar.BEC_transformed_equivalent(epsilon, idx_remain, N-1);
if idx_current == 0
    % bad channel
    y = 2 * y - y^2;
else
    % good channel
    y = y^2;
end


end