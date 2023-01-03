% 4QAM 

function modulated_signal = Sixteen_QAM(coded_signal, Ebav)
    switch_int = 0;
    for i = 1 : 4
        coded_signal(i)
        switch_int = switch_int + 2^(4-i)*coded_signal(i)
    switch switch_int
        case 0
            modulated_signal = [sqrt(Ebav/5) sqrt(Ebav/5)];
        case 1
            modulated_signal = [sqrt(Ebav/5) 2*sqrt(Ebav/5)];
        case 2
            modulated_signal = [2*sqrt(Ebav/5) sqrt(Ebav/5)];
        case 3
            modulated_signal = [2*sqrt(Ebav/5) 2*sqrt(Ebav/5)];
        case 4
            modulated_signal = [sqrt(Ebav/5) -sqrt(Ebav/5)];
        case 5
            modulated_signal = [sqrt(Ebav/5) -2*sqrt(Ebav/5)];
        case 6
            modulated_signal = [2*sqrt(Ebav/5) -sqrt(Ebav/5)];
        case 7
            modulated_signal = [2*sqrt(Ebav/5) -2*sqrt(Ebav/5)];
        case 8
            modulated_signal = [-sqrt(Ebav/5) sqrt(Ebav/5)];
        case 9
            modulated_signal = [-sqrt(Ebav/5) 2*sqrt(Ebav/5)];
        case 10
            modulated_signal = [-2*sqrt(Ebav/5) sqrt(Ebav/5)];
        case 11
            modulated_signal = [-2*sqrt(Ebav/5) 2*sqrt(Ebav/5)];
        case 12
            modulated_signal = [-sqrt(Ebav/5) -sqrt(Ebav/5)];
        case 13
            modulated_signal = [-sqrt(Ebav/5) -2*sqrt(Ebav/5)];
        case 14
            modulated_signal = [-2*sqrt(Ebav/5) -sqrt(Ebav/5)];
        case 15
            modulated_signal = [-2*sqrt(Ebav/5) -2*sqrt(Ebav/5)];
    end
end