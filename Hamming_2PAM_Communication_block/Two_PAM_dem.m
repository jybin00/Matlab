% 2-PAM demodulator LLR(Log Likelihood Ratio)


function demoulated_output = Two_PAM_dem(received_signal)     % 채널을 통과한 신호와 노이즈의 variance를 받는다.
    for i = 1:length(received_signal)
        if received_signal(i) > 0
            received_signal(i) = 1;
        else                                     % 0보다 작으면 0으로 판단
            received_signal(i) = 0;
        end
    end
    demoulated_output = received_signal;
    
end