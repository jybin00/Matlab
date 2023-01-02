%FER analysis

function correct_or_error = FER_analysis(input, estimation)
    correct_or_error = [0 0];
    if input == estimation
        correct_or_error(1) = correct_or_error(1) + 1;
    else
        correct_or_error(2) = correct_or_error(2) + 1;
    end
end