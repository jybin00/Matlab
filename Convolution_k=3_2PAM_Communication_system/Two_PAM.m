%% 2PAM (Pulse Amplitude Modulation) Modulator

function modulated_signal = Two_PAM(code, Es)   % Original code, bit energy
%     for i = 1:length(code)
%         if code(i) == 1
%             code(i) = sqrt(Es);
%         else
%             code(i) = -sqrt(Es);
%         end
%     end
%     modulated_signal = code;
modulated_signal = sqrt(Es) * (2*code -1);
end