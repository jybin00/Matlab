function y = binary_erasure_channel(p, x)
flip = rand(size(x)) < p;
y = x;
y(flip) = nan;
end