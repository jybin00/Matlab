function output_sym = QPSK_Modulation(N_sym,input_bit) 

for i=1:N_sym
    if input_bit(2*i-1)==1 && input_bit(2*i)==1
        output_sym(i) = (1/sqrt(2))*(1+1i);
    elseif input_bit(2*i-1)==1 && input_bit(2*i)==0
        output_sym(i) = (1/sqrt(2))*(1-1i);
    elseif input_bit(2*i-1)==0 && input_bit(2*i)==1
        output_sym(i) = (1/sqrt(2))*(-1+1i);
    elseif input_bit(2*i-1)==0 && input_bit(2*i)==0
        output_sym(i) = (1/sqrt(2))*(-1-1i);
    end 
end
end