

function output = bin_deci (input)
    num = length(input);
    d=0;
    for i = 1:  num
        d = d + (2^(num-i))*input(1,i);
    end
    output = d;
end