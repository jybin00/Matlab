% ML demodulation
% ML Metric이 가장 작은 값으로 demodulation

function demodulated_output = Demodulation(received_signal)
%     if received_signal > 0                     
%         demodulated_output = 1;
%     else
%         demodulated_output = 0;
%     end
demodulated_output  = received_signal > 0;
end
