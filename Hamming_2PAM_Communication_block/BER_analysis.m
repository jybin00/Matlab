% BER 분석하기.

function correct_or_error = BER_analysis(input, estimation)
    correct_or_error = [0 0];
    for i = 1:length(estimation)
        if input(i) == estimation(i) 
            correct_or_error(1) = correct_or_error(1) + 1;
        else
            correct_or_error(2) = correct_or_error(2) + 1;
        end
    end
end