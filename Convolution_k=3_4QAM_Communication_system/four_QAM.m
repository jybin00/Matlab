% 4QAM 
% Grey coding 적용
% dmin = 2*sqrt(Ebav/2)

function modulated_signal = four_QAM(coded_signal, Es)
    switch_int = 0;
    for i = 1 : 2
        switch_int = switch_int + 2^(2-i)*coded_signal(i);
        switch switch_int               % grey coding
            case 0                          % [0 0]
                modulated_signal = sqrt(Es/2)+ 1i*sqrt(Es/2);
            case 1                          % [0 1]
                modulated_signal = -sqrt(Es/2)+1i* sqrt(Es/2);
            case 2                          % [1 0]
                modulated_signal = sqrt(Es/2)-1i*sqrt(Es/2);
            case 3                          % [1 1]
                modulated_signal = -sqrt(Es/2)-1i*sqrt(Es/2);
        end
    end
end