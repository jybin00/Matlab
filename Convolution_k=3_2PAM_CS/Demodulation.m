% ML demodulation
% ML Metric이 가장 작은 값으로 demodulation

function demodulated_output = Demodulation(received_signal)
    if real(received_signal) > 0                                   % 우반면
        if imag(received_signal) >0                                 % 1사분면   % 00
            demodulated_output = [0 0];
        else                                                             % 4사분면
            demodulated_output = [1 0];
        end
    else                                                                % 좌반면
        if imag(received_signal) > 0                             % 2사분면
            demodulated_output = [0 1];
        else
            demodulated_output = [1 1];                 % 3사분면
        end
    end
end