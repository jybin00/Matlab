function output_bit = Demodulation(M,N_sym,input_symbol)

for k=1:N_sym
    if M==1
        if input_symbol(k) < 0
            output_bit(k) = 0;
        elseif input_symbol(k) >= 0
            output_bit(k) = 1;
        end
    elseif M==2
        if real(input_symbol(k)) >= 0 && imag(input_symbol(k)) >= 0 
            output_bit(2*k-1) = 0;
            output_bit(2*k) = 0;
    elseif real(input_symbol(k)) < 0 && imag(input_symbol(k)) >= 0 
        output_bit(2*k-1) = 1;
        output_bit(2*k) = 0;
    elseif real(input_symbol(k)) < 0 && imag(input_symbol(k)) < 0 
        output_bit(2*k-1) = 1;
        output_bit(2*k) = 1;
    elseif real(input_symbol(k)) >= 0 && imag(input_symbol(k)) < 0
            output_bit(2*k-1) = 0;
            output_bit(2*k) = 1;
        end 
    end
end 
end