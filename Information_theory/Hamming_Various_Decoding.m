% Hamming code % Peeling decoder
% FER measure at BEC according to SNR
% error rate 0~0.5
clc
tic

% SNR & Noise variance setting
% SNR_dB = 0:10;
% N0 = 1./(10.^(SNR_dB./10)

% Number of message "frame"
num_f = 100000;

% Error rate
error_r = 0.0:0.01:1;

% FER var setting || p = peeling, s = syndrome, M = ML
FER_p = zeros(1, length(error_r));
FER_s = zeros(1, length(error_r));
FER_M = zeros(1, length(error_r));
FER_M_BEC = zeros(1, length(error_r));

% message generation
message = randi([0, 1], [num_f, 4]);

% G matrix
G = [1 1 1 0 0 0 0;
       1 0 0 1 1 0 0;
       0 1 0 1 0 1 0;
       1 1 0 1 0 0 1;];

% for loop
for j = 1:length(error_r)
    fprintf("error rate : %d\n", (error_r(j)))
    for i = 1: num_f
        m = message(i, :);
        Coded_frame = mod(G' * message(i, : )', 2);
        Coded_frame_view = Coded_frame';

        err = err_frame(Coded_frame', error_r(j));

        de_m1 = Peeling(err(1, :));
        de_m2 = Syndrome(err(2, :));
        de_m3 = ML_decoding(err(3, :));
        de_m4 = ML_decoding_BEC(err(1, :));

        is_err = error_check (message(i, :)', de_m1, de_m2, de_m3, de_m4);
        FER_p(j) = FER_p(j) + is_err(1);
        FER_s(j) = FER_s(j)+ is_err(2);
        FER_M(j) = FER_M(j) + is_err(3);
        FER_M_BEC(j) = FER_M_BEC(j) + is_err(4);
    end
end
toc
% decoding

%%
% plot
close all
figure(1);
hold on
grid on
set(gca, "YScale", "log")
plot(error_r, FER_s/ num_f, "Linewidth", 2, "Marker", "+")
plot(error_r, FER_M/ num_f, "Linewidth", 2,  "LineStyle", "--", "Marker", "o")
xlabel('epsilon')
ylabel('FER')
legend("Syndrome decoding", "ML decoding")

figure(2);
hold on
grid on
set(gca, "YScale", "log")
plot(error_r, FER_p/ num_f, "Linewidth", 2, "Marker", "+")
plot(error_r, FER_M_BEC/ num_f, "Linewidth", 2, "Marker", "o")
xlabel('epsilon')
ylabel('FER')
legend("Peeling decoder", "BEC ML decoder")

%%

% Peeling decoder
function decoding = Peeling(Frame)
    while(true)
        Frame_past = Frame;

        a = Frame([1,3,5,7]);
        if sum(a == 2) == 1
            if mod(sum(a(a~=2)), 2) == 0
                a(a == 2) = 0;
            else
                a(a == 2) = 1;
            end
            Frame([1,3,5,7]) = a;
        end

        b = Frame([2,3,6,7]);
        if sum(b == 2) == 1
            if mod(sum(b(b~=2)), 2) == 0
                b(b == 2) = 0;
            else
                b(b == 2) = 1;
            end
            Frame([2,3,6,7]) = b;
        end
    
        c = Frame([4,5,6,7]);
        if sum(c == 2) == 1
            if mod(sum(c(c~=2)), 2) == 0
                c(c == 2) = 0;
            else
                c(c == 2) = 1;
            end
            Frame([4,5,6,7]) = c;
        end
        if(Frame_past == Frame)
            break
        end
    end
    decoding = Frame([3,5,6,7]);
end

% Syndrome decoder
function decoding = Syndrome(Frame)
H = [0 0 0 1 1 1 1;
    0 1 1 0 0 1 1;
    1 0 1 0 1 0 1];
e = mod(H*Frame', 2);
error_position = e(1)*2^2 + e(2)*2^1 + e(3)*2^0;
if error_position > 0
    e_c = zeros(1,7);
    e_c(error_position) = 1;
    error_correction = mod(Frame + e_c, 2);
    decoding = error_correction([3,5,6,7]);
else
    decoding = Frame([3, 5, 6, 7]);
end
end

% ML decoder
function decoding = ML_decoding(Frame)

code_word = [0 0	0	0	0	0	0;
    1	1	0	1	0	0	1;
    0	1	0	1	0	1	0;
    1	0	0	0	0	1	1;
    1	0	0	1	1	0	0;
    0	1	0	0	1	0	1;
    1	1	0	0	1	1	0;
    0	0	0	1	1	1	1;
    1	1	1	0	0	0	0;
    0	0	1	1	0	0	1;
    1	0	1	1	0	1	0;
    0	1	1	0	0	1	1;
    0	1	1	1	1	0	0;
    1	0	1	0	1	0	1;
    0	0	1	0	1	1	0;
    1	1	1	1	1	1	1];
ML_metric = zeros(16, 1);
Error_Sqr = (code_word - Frame).^2;  % 이렇게 계산하면 모든 행에서 값이 빠짐.
for j = 1:16
    ML_metric(j, 1) = sum(Error_Sqr(j, :));                   % ML_metric 계산. 코드에서 받은 값 빼서 제곱.
end
[~, min_distance] = min(ML_metric);     % 최소가 되는 index 찾기.
decoding_output = code_word(min_distance, :);   % 찾은 index로 decoding하기.
decoding = decoding_output([3,5,6,7]);
end

% Error Frame generator.
function err_frame = err_frame(Coded_frame, error_rate)

err_frame = zeros(3,7);

a = rand([1, 7]) < error_rate;
err_frame(1, :) =  Coded_frame;
err_frame(2, :) =  Coded_frame;
err_frame(3, :) =  Coded_frame;

err_frame(2, a == 1) = ~Coded_frame(a == 1);
err_frame(3, a == 1) = ~Coded_frame(a == 1);
Coded_frame(a == 1) = 2;
err_frame(1,:) = Coded_frame;

end

% Different decoded message makes an error.
function is_err = error_check(msg, decoded_msg1, decoded_msg2, decoded_msg3, decoded_msg4)

msg = msg';
is_err = zeros(1,4);
    if(msg == decoded_msg1)
        is_err(1) = 0;
    else
        is_err(1) = 1;
    end
    
    if(msg == decoded_msg2)
        is_err(2) = 0;
    else
        is_err(2) = 1;
    end
    
    if(msg == decoded_msg3)
        is_err(3) = 0;
    else
        is_err(3) = 1;
    end

    if(msg == decoded_msg4)
        is_err(4) = 0;
    else
        is_err(4) = 1;
    end
end

% ML decoder
function decoding = ML_decoding_BEC(Frame)

code_word = [0 0	0	0	0	0	0;
    1	1	0	1	0	0	1;
    0	1	0	1	0	1	0;
    1	0	0	0	0	1	1;
    1	0	0	1	1	0	0;
    0	1	0	0	1	0	1;
    1	1	0	0	1	1	0;
    0	0	0	1	1	1	1;
    1	1	1	0	0	0	0;
    0	0	1	1	0	0	1;
    1	0	1	1	0	1	0;
    0	1	1	0	0	1	1;
    0	1	1	1	1	0	0;
    1	0	1	0	1	0	1;
    0	0	1	0	1	1	0;
    1	1	1	1	1	1	1];

a = (Frame ~= 2);
ML_metric = zeros(16, 1);
Error_Sqr = (code_word(:, a) - Frame(a)).^2;  % 이렇게 계산하면 모든 행에서 값이 빠짐.
for j = 1:16
    ML_metric(j, 1) = sum(Error_Sqr(j, :));                   % ML_metric 계산. 코드에서 받은 값 빼서 제곱.
end
[~, min_distance] = min(ML_metric);     % 최소가 되는 index 찾기.
decoding_output = code_word(min_distance, :);   % 찾은 index로 decoding하기.
decoding = decoding_output([3,5,6,7]);
end









