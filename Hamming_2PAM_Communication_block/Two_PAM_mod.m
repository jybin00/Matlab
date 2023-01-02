%% 2PAM (Pulse Amplitude Modulation) Modulator

function modulated_signal = Two_PAM_mod(code, Eb)   % Original code, bit energy
    for i = 1:length(code)
        if code(i) == 1
            code(i) = sqrt(Eb);
        else
            code(i) = -sqrt(Eb);
        end
    end
    modulated_signal = code;
end