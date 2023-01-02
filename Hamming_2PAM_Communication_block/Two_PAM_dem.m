% 2-PAM demodulator LLR(Log Likelihood Ratio)


function demoulated_output = Two_PAM_dem(received_signal, sigma_v, Eb)     % 채널을 통과한 신호와 노이즈의 variance를 받는다.
    for i = 1:length(received_signal)
        l_po= (1/(sigma_v*sqrt(2*pi)))*exp(-((received_signal(i)-sqrt(Eb))^2)/(2*sigma_v^2));    % positive one
        l_no =  (1/(sigma_v*sqrt(2*pi)))*exp(-((received_signal(i)+sqrt(Eb))^2)/(2*sigma_v^2));  % negative one
    
        if log(l_po/l_no) > 0   % ML로 가정하고 LLR 계산해서 0보다 크면 1로 판단.
            received_signal(i) = 1;
        else                            % 0보다 작으면 0으로 판단
            received_signal(i) = 0;
        end
    end
    demoulated_output = received_signal;
    
end