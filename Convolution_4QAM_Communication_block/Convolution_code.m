% Convolution coding 

function coded_signal = Convolution_code(input)
    m = zeros(6, 1);
    for i = 1:length(input)
        coded_signal(2*i-1) = mod(input(i) + m(1) + m(2) + m(3),2);
        coded_signal(2*i) = mod(input(i) + m(1) + m(2) + m(3) + m(5) + m(6), 2);
        m(1) = input(i);
        for j = 2 : 6
            m(j) = m(j-1);
        end
    end
end