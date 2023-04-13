function output_sym = Modulation(M,N_sym,input_bit)

for i=1:N_sym
    if M==1
        if input_bit(i)==0
            output_sym(i) = -1;
        elseif input_bit(i)==1
            output_sym(i) = 1;
        end
    elseif M==2
        if input_bit(2*i-1)==0 && input_bit(2*i)==0 
            output_sym(i) = (1/sqrt(2))*(1+1i);
        elseif input_bit(2*i-1)==0 && input_bit(2*i)==1 
            output_sym(i) = (1/sqrt(2))*(1-1i);
        elseif input_bit(2*i-1)==1 && input_bit(2*i)==0 
            output_sym(i) = (1/sqrt(2))*(-1+1i);
        elseif input_bit(2*i-1)==1 && input_bit(2*i)==1 
            output_sym(i) = (1/sqrt(2))*(-1-1i);
        end 
    end
end 
end