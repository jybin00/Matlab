% BER 분석하기.

function correct_or_error = BER_analysis(input, estimation)  % compare with input and output
    correct_or_error = [0 0];                                                   % correct or error array
    for i = 1:length(estimation)
        if input(i) == estimation(i)                                            % if input is same with estimation
            correct_or_error(1) = correct_or_error(1) + 1;          % increase correct index
        else
            correct_or_error(2) = correct_or_error(2) + 1;           % else increase error index
        end
    end
end