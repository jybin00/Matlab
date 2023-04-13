function output_sym = Rayleigh_Fading(N_sym,input_sym)

h_Rayleigh = abs(1/sqrt(2)*(randn(1, N_sym) +1i*randn(1, N_sym)));
output_sym = h_Rayleigh.*input_sym;

end