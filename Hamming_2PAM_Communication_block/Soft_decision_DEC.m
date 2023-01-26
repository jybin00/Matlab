%% Haming code soft decision decoding

function decoding_output = Soft_decision_DEC(received_signal, Es)
    modified_code = [-1  -1  -1  -1  -1   -1  -1;
                                  1	1	-1	1	-1  -1	  1;
                                -1	1	-1	1	-1	1	-1;
                                  1  -1	-1  -1	-1	1	  1;
                                  1  -1	-1	1	  1  -1	-1;
                                -1	1	-1  -1	  1  -1	  1;
                                  1	1	-1  -1	  1	1	-1;
                                -1  -1	-1	1	  1	1	  1;
                                  1	1	  1  -1	-1  -1	-1;
                                -1  -1	  1	1	-1  -1	  1;
                                  1  -1	  1	1	-1	1	-1;
                                -1	1	  1  -1	-1	1	  1;
                                -1	1	  1	1	  1  -1	-1;
                                  1  -1	  1  -1	  1  -1	  1;
                                -1  -1	  1  -1	  1	1	-1;
                                  1	1	  1	1	  1	1	  1];

    modified_code = modified_code * sqrt(Es);

    code_word = [0	0	0	0	0	0	0;
                        1	1	0	1	0	0	1;
                        0	1	0	1	0	1	0;
                        1	0	0	0	0	1	1;
                        1	0	0	1	1	0	0;
                        0	1	0	0	1	0	1;
                        1	1	0	0	1	1	0;
                        0	0	0	1	1	1	1;
                        1	1	1	0	0	0	0;
                        0	0	1	1	0	0	1;
                        1	0	1	1	0	1	0;
                        0	1	1	0	0	1	1;
                        0	1	1	1	1	0	0;
                        1	0	1	0	1	0	1;
                        0	0	1	0	1	1	0;
                        1	1	1	1	1	1	1];
    ML_metric = zeros(16, 1);
    Error_Sqr = (modified_code - received_signal).^2;  % 이렇게 계산하면 모든 행에서 값이 빠짐.
    for j = 1:16
        ML_metric(j, 1) = sum(Error_Sqr(j, :));                   % ML_metric 계산. 코드에서 받은 값 빼서 제곱.
    end
    [~, min_distance] = min(ML_metric);     % 최소가 되는 index 찾기.
    decoding_output = code_word(min_distance, [3, 5, 6, 7]);   % 찾은 index로 decoding하기.
end