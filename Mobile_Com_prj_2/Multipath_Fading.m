function output_sym = Multipath_Fading(N_sym,input_sym,h_tap) 

output_sym = cconv(input_sym, h_tap, N_sym);

end