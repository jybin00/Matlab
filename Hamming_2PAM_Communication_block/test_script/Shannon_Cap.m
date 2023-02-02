%C = argmax(px) (H(y) - H(y|x)) function 
% px = 1/2 가정 

function Capacity = Shannon_Cap (eb, sigma_v)
    syms x
    coefficient = 1/(sqrt(8*pi*sigma_v^2));
    exponential = (exp(-(x-sqrt(abs(eb))).^2/(2*sigma_v^2))+exp(-(x+sqrt(abs(eb))).^2/(2*sigma_v^2)));
    fy = coefficient * exponential;
    minus_fylogfy = -1*fy*log2(fy);

    Capacity = vpaintegral(minus_fylogfy, x, [-inf inf]) - 1/2*(log2(2*pi*exp(1)*sigma_v^2));
end