%ML demodulation

function demodulated_output = ML_demodulation(received_signal, Esav)
    if received_signal(1) > 0                                   % 우반면
        if received_signal(2)>0                                 % 1사분면   % 0001    0011
            if received_signal(1) > 2*sqrt(Esav/10)                     % 0000    0010
                if received_signal(2) > 2*sqrt(Esav/10)  
                    demodulated_output = [0 0 1 1];
                else
                    demodulated_output = [0 0 1 0];
                end
            else
                if received_signal(2) > 2*sqrt(Esav/10)
                    demodulated_output = [0 0 0 1];
                else
                    demodulated_output = [0 0 0 0];
                end
            end
        else                                                             % 4사분면   % 0100   0110
            if received_signal(1) > 2*sqrt(Esav/10)                    % 0101   0111
                if received_signal(2) < -2*sqrt(Esav/10)
                    demodulated_output = [0 1 1 1];
                else
                    demodulated_output = [0 1 1 0];
                end
            else
                if received_signal(2) < -2*sqrt(Esav/10)
                    demodulated_output = [0 1 0 1];
                else
                    demodulated_output = [0 1 0 0];
                end
            end
        end
    else                                                                % 좌반면
        if received_signal(2) > 0                             % 2사분면   % 1011   1001
            if received_signal(1) < -2*sqrt(Esav/10)                 %  1010  1000
                if received_signal(2) > 2*sqrt(Esav/10)
                    demodulated_output = [1 0 1 1];
                else
                    demodulated_output = [1 0 1 0];
                end
            else
                if received_signal(2) > 2*sqrt(Esav/10)
                    demodulated_output = [1 0 0 1];
                else
                    demodulated_output = [1 0 0 0];
                end
            end
        else                                                        % 3사분면    % 1110   1100
            if received_signal(1) < -2*sqrt(Esav/10)               % 1111   1101
                if received_signal(2) < -2*sqrt(Esav/10)
                    demodulated_output = [1 1 1 1];
                else
                    demodulated_output = [1 1 1 0];
                end
            else
                if received_signal(2) < -2*sqrt(Esav/10)
                    demodulated_output = [1 1 0 1];
                else
                    demodulated_output = [1 1 0 0];
                end
            end
        end
    end
end