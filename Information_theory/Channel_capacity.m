clear
close all
syms x y 

SNR_dB = -10:30;
SNR_lin = 10.^(SNR_dB./10);
sigma = sqrt(1./SNR_lin);
h_y = zeros(length(SNR_dB), 1);
I_x_y = zeros(length(SNR_dB), 2);
h_y_x = zeros(length(SNR_dB), 2);
% mod 1 => BPSK, 2 => 4PAM, 3 => both
mod = 4;


% p(y) = p(x=1)*p(y|x=1) + p(x=-1)*p(y|x=-1)
if mod == 1||3||4
    p_p1 = 1/2;
    p_m1 = 1/2;
    for i = 1:length(SNR_dB) 
        p(x,y) = 1./sqrt(2*pi*sigma(i).^2).*exp(-(y-x)^2./(2*sigma(i).^2));
        p_y = p_m1*p(y,-1) + p_p1*p(y, 1);
        h_y(i) = vpaintegral(p_y*log2(1/p_y), y, [-inf, inf]);
        h_y_x(i,1) = 1/2*log2(2*pi*exp(1)*sigma(i)^2);
        I_x_y(i,1) = h_y(i) - h_y_x(i, 1);
    end
end

if mod == 2||3
    p_m3 = 1/4;
    p_m1 = 1/4;
    p_p1 = 1/4;
    p_p3 = 1/4;
    sigma = sqrt(5./SNR_lin);
    for i = 1:length(SNR_dB) 
        p(x,y) = 1./sqrt(2*pi*sigma(i).^2).*exp(-(y-x)^2./(2*sigma(i).^2));
        p_y = p_m3*p(y,-3) + p_m1*p(y,-1) + p_p1*p(y,1) + p_p3*p(y,3);
        h_y(i) = vpaintegral(p_y*log2(1/p_y), y, [-inf, inf]);
        h_y_x(i,2) = 1/2*log2(2*pi*exp(1)*sigma(i)^2);
        I_x_y(i,2) = h_y(i) - h_y_x(i, 2);
    end
end
SNR_dB;
x = 1/2*log2(1+SNR_lin);
figure

%%
close
plot(SNR_dB, x, 'r')
hold on

if mod == 1
    plot(SNR_dB, I_x_y(:, 1), 'b')
    legend('Shannon capacity limit', 'BPSK')
end

if mod == 2
    plot(SNR_dB, I_x_y(:, 2), 'b')
    legend('Shannon capacity limit', '4PAM')
end

if mod == 3
    plot(SNR_dB, I_x_y(:, 2), 'b')
    hold on
    plot(SNR_dB, I_x_y(:, 1), 'g')
    legend('Gaussian input', '4PAM', 'BPSK')
end

if mod == 4
    plot(SNR_dB, I_x_y(:, 1), 'b')
    hold on
    pe = qfunc(sqrt(SNR_lin));
    pep = 1-pe;
    H_pe = -pe.*log2(pe) - (pep).*log2(pep);
    plot(SNR_dB, 1- H_pe, 'g');
    legend('$\frac{1}{2} log_2(1+SNR)$', 'BPSK', '$1-H(Q(\sqrt{SNR})$', 'Interpreter', 'latex')
end
axis([-10 30 0 6])
title('I(X;Y) (bit/channel use) vs SNR (dB)')
xlabel('SNR (dB)');
ylabel('I(X;Y) (bit/channel use)')





