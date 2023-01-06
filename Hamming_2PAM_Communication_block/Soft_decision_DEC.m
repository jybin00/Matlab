%% Haming code soft decision decoding

function decoded_ouptut = Soft_decision_DEC(received_signal)
    code = [0	0	0	0	0	0	0;
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
    Error_Sqr = (code - received_signal).^2;  % 이렇게 계산하면 모든 행에서 값이 빠짐.
    for i = 1:16
        for j = 1:7
            ML_metric(i, 1) = ML_metric(i, 1) + Error_Sqr(i, j);   % ML_metric 계산. 코드에서 받은 값 빼서 제곱.
        end
    end
    [~, min_distance] = min(ML_metric);     % 최소가 되는 index 찾기.
    decoded_ouptut = code(min_distance, [3 5 6 7]);   % 찾은 index로 decoding하기.
end