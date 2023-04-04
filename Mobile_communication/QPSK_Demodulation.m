function output_bit = QPSK_Demodulation(N_sym,input_symbol)

for k=1:N_sym
    if real(input_symbol(k)) >= 0 && imag(input_symbol(k)) >= 0
        output_bit(2*k-1) = 1;
        output_bit(2*k) = 1;
    elseif real(input_symbol(k)) < 0 && imag(input_symbol(k)) >= 0
        output_bit(2*k-1) = 0;
        output_bit(2*k) = 1;
    elseif real(input_symbol(k)) < 0 && imag(input_symbol(k)) < 0
        output_bit(2*k-1) = 0;
        output_bit(2*k) = 0;
    elseif real(input_symbol(k)) >= 0 && imag(input_symbol(k)) < 0
        output_bit(2*k-1) = 1;
        output_bit(2*k) = 0;
    end
end
