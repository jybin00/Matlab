function G = generator_matrix(n)
% ========================================
% Polar code를 위한 generator matrix 반환
% Input:
%   n: matrix size (power of 2...)
% Output:
%   G: generator matrix
% ========================================

assert (n > 0)
assert (abs(log2(n) - floor(log2(n))) < 1e-10)


G = 1;
F = [1, 0; 1, 1];
idx = log2(n);

for i = 1:idx
    G = kron(F, G);
end

% G = G((bitrevorder(0:n-1)+1),:);

end