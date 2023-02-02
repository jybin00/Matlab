% Convolution coding (3, [8, 5]) = (3, [111], [101])

function coded_signal = Convolution_code(input)     % input sequence를 받는다.
    coded_signal = zeros(1, (length(input)+2)*2);
    m = zeros(2, 1);                                                    % memory를 0으로 초기화 한다.
    input = [input zeros(1, 2)];                                    % message bit 후에 memory initialization을 위해서 zero padding
    
    for i = 1:length(input)                                                                             % input message의 길이만큼 반복문 실행
        coded_signal(1, 2*i-1, 1) = mod(input(i) + m(1) + m(2),2);      % 1+1+1
        coded_signal(1,2*i) = mod(input(i) + m(2), 2);                         % 1+0+1
        m(2) = m(1);                                                                            % register shift
        m(1) = input(i);                                                                        % input m1으로 이동
    end
end