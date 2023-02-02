% 4QAM 

function modulated_signal = Sixteen_QAM(coded_signal, Esav)
    switch_int = 0;
    for i = 1 : 4
        coded_signal(i);
        switch_int = switch_int + 2^(4-i)*coded_signal(i);
        switch switch_int               % grey coding
            case 0                          % [0 0 0 0]
                modulated_signal = sqrt(Esav/10) +1i*sqrt(Esav/10);
            case 1                          % [0 0 0 1]
                modulated_signal = sqrt(Esav/10) +1i*3*sqrt(Esav/10);
            case 2                          % [0 0 1 0]
                modulated_signal = 3*sqrt(Esav/10) +1i*sqrt(Esav/10);
            case 3                          % [0 0 1 1]
                modulated_signal = 3*sqrt(Esav/10) +1i*3*sqrt(Esav/10);
            case 4                          % [0 1 0 0]
                modulated_signal = sqrt(Esav/10) -1i*sqrt(Esav/10);
            case 5                          % [0 1 0 1]
                modulated_signal = sqrt(Esav/10) -1i*3*sqrt(Esav/10);
            case 6                          % [0 1 1 0]
                modulated_signal = 3*sqrt(Esav/10) -1i*sqrt(Esav/10);
            case 7                          % [0 1 1 1]
                modulated_signal = 3*sqrt(Esav/10) -1i*3*sqrt(Esav/10);
            case 8                          % [1 0 0 0]
                modulated_signal = -sqrt(Esav/10) +1i*sqrt(Esav/10);
            case 9                          % [1 0 0 1]
                modulated_signal = -sqrt(Esav/10) +1i*3*sqrt(Esav/10);
            case 10                        % [1 0 1 0]
                modulated_signal = -3*sqrt(Esav/10) +1i*sqrt(Esav/10);
            case 11                         % [1 0 1 1]
                modulated_signal = -3*sqrt(Esav/10) +1i*3*sqrt(Esav/10);
            case 12                         % [1 1 0 0]
                modulated_signal = -sqrt(Esav/10) -1i*sqrt(Esav/10);
            case 13                         % [1 1 0 1]
                modulated_signal = -sqrt(Esav/10) -1i*3*sqrt(Esav/10);
            case 14                         % [1 1 1 0]
                modulated_signal = -3*sqrt(Esav/10) -1i*sqrt(Esav/10);
            case 15                         % [1 1 1 1]
                modulated_signal = -3*sqrt(Esav/10) -1i*3*sqrt(Esav/10);
        end
    end
end