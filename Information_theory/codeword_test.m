
code_word = [0 0	0	0	0	0	0;
    1	1	0	1	0	0	1;
    0	1	0	1	0	1	0;
    1	0	0	0	0	1	1;
    1	0	0	1	1	0	0;
    0	1	0	0	1	0	1;
    1	1	0	0	1	1	0;
    0	0	0	1	1	1	1;
    1	1	1	0	0	0	0;
    0	0	1	1	0	0	1;
    1	0	1	1	0	1	0;
    0	1	1	0	0	1	1;
    0	1	1	1	1	0	0;
    1	0	1	0	1	0	1;
    0	0	1	0	1	1	0;
    1	1	1	1	1	1	1];

i = zeros(1,16);
G = [1 1 1 0 0 0 0;
    1 0 0 1 1 0 0;
    0 1 0 1 0 1 0;
    1 1 0 1 0 0 1;];

code = zeros(16, 7);

for j = 1:16
    u = dec2bools(j-1, 4);
    code(j, :) = mod(u*G, 2);
    
end
a = code_word - code;

function bools = dec2bools(dec,len)
    bools = zeros([1,len]);
    bools_nopad = (dec2bin(dec)-'0');
    nopad_len = length(bools_nopad);
    bools(len-nopad_len+1:end) = bools_nopad;
end