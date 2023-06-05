function y = AWGN_transformed_equivalent(p, idx, N)
% https://publik.tuwien.ac.at/files/publik_262980.pdf 
% Algorithm 3
% p: qfunc(snr), snr in linear scale (not dB)

idx_bin = de2bi(idx, N);

assert (~isempty(idx_bin))
assert (N > 0)


idx_current = idx_bin(1);
if length(idx_bin) == 1
    if idx_current == 0
        % bad channel
        y = 2*p*(1-p);
    else
        % good channel
        y = qfunc(sqrt(2)*qfuncinv(p));
    end
    return
end

idx_remain = bi2de(idx_bin(2:end));
y = polar.AWGN_transformed_equivalent(p, idx_remain, N-1);
if idx_current == 0
    % bad channel 
    y = 2*y*(1-y);
else
    % good channel
    y = qfunc(sqrt(2)*qfuncinv(y));
end


end