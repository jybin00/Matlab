% Convolution coding (7, [171, 133]) = (7, [1111001], [1011011])

function coded_signal = Convolution_code_tail_bit(input)     % input sequence를 받는다.
    nmem = 6;
    coded_signal = zeros(1, length(2*(input+nmem)));
    m = zeros(nmem, 1);                                                    % memory를 0으로 초기화 한다.
    input = [input zeros(1, nmem)];                                    % message bit 후에 memory initialization을 위해서 zero padding
    
    for i = 1:length(input)                                                                             % input message의 길이만큼 반복문 실행
        coded_signal(1, 2*i-1, 1) = mod(input(i) + m(1) + m(2) + m(3) + m(6),2);             % 1+1+1+1+0+0+1
        coded_signal(1,2*i) = mod(input(i) + m(2) + m(3) + m(5) + m(6), 2);     % 1+0+1+1+0+1+1
        for j = 2 : 6
            m(8-j) = m(7-j);                                                                              % register shift
        end
        m(1) = input(i);                                                                                    % input m1으로 이동
    end
end