clc
clear all
close all

%% Simulation in multipath channel for delay that exceeds the guard time by 100 percents of the fft interval
L_d = 200; % Length of delay
L_CP = 100; % Length of CP


%% System parameter setting, memory allocation
NumSymbol = 1024;
SNR_dB = 50;

h = zeros(1, NumSymbol);
xn_cp = zeros(1, NumSymbol + L_CP);
y = zeros(1, NumSymbol);
Xk_det = zeros(1, NumSymbol);


%% Transmitter
% Generate random signal
Xk = (1/sqrt(2))*(2*(randi(2, [1, NumSymbol]) - 1.5) + 1j * 2*(randi(2, [1, NumSymbol]) - 1.5));

% Perform IFFT to get time-domain signal
xn = ifft(Xk);

% Add cyclic prefix
for a = 1:L_CP
    xn_cp(a) = xn(NumSymbol - L_CP + a);
end

for b = (L_CP+1):NumSymbol+L_CP
    xn_cp(b) = xn(b - L_CP);
end


%% Multipath Channel
% Generate multipath channel with delay
h = [linspace(1, 0.75, (L_d+1)), zeros(1, NumSymbol - (L_d+1))];

% Generate Rayleigh fading channel
h_fading = (randn(1, NumSymbol) + 1j*randn(1, NumSymbol))/sqrt(2);

% Apply channel to signal with fading
yn = conv(xn_cp, h.*h_fading);


%% Receiver
% Remove cyclic prefix
for k = 1:NumSymbol
    y(k) = yn(k+L_CP);
end

% Add 25 dB noise to the received signal
noise_power = 10^(-SNR_dB/10); % Convert noise power to linear scale noise = sqrt(noise_power/2) * (randn(1, NumSymbol) + 1j*randn(1, NumSymbol)); % Generate complex Gaussian noise
y = y + noise; % Add noise to received signal

% Perform FFT on the received signal
Yk = fft(y);

% Calculate channel frequency response
Hk = fft(h.*h_fading);

% Perform subcarrier-wise equalization
Xk_det = Yk./Hk; %  Implement properly! 

% Calculate channel gain in dB
channel_gain_dB = 10*log10(abs(Hk));

% Define frequency axis in fraction of sampling frequency
f = linspace(-0.5, 0.5, NumSymbol);


%% Simulation Result
% Plot received constellation diagram
scatter(real(Xk_det), imag(Xk_det))
axis([-3.5 3.5 -3.5 3.5])
grid on
box on
title(['OFDM simulation with Rayleigh fading, channel delay exceeds guard time by ', ...
    num2str((L_d/L_CP)*100), '%'],'interpreter','latex')


% Plot channel gain in frequency domain
figure
plot(f, channel_gain_dB)
grid on
box on
title(['Channel gain in frequency domain with Rayleigh fading: channel delay ', ...
    num2str(L_d),'*Ts'],'interpreter','latex')
xlabel('Normalized frequency')
ylabel('Magnitude (dB)')