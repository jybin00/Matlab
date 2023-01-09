% Convolution coding (3, [111, 101])

function coded_signal = Convolution_code_tail_bit_X(input, memory)     % input sequence를 받는다. % memory 값을 받는다.
    
    for i = 1:length(input)                                                                             % input message의 길이만큼 반복문 실행
        coded_signal(2*i-1) = mod(input(i) + memory(1) + memory(2) ,2);             % 1 + 1 + 1
        coded_signal(2*i) = mod(input(i) + memory(2), 2);     % 1 + 0 + 1
        memory(2) = memory(1);                                                                              % register shift
        memory(1) = input(i);                                                                          % input m1으로 이동
    end
end