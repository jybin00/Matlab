syms x
idf = zeros(length(1:0.1:5.4),1);

H = @(p) -p*log2(p) -(1-p)*log2(1-p);
k = 1;
i = 5.3;
C= Shannon_Cap(i, 2)
a= (4-7*a)/4
p = H(10^-7)
idf = solve(Shannon_Cap(x, 2) == p, x)
real(idf)
10*log10(i/8)
%%
close
i = 1:0.1:5.4;
db = 10*log10(i/8);

scatter(db, 10.^log10(idf));
hold on
grid,axis([-2 14 10^(-7) 10^-(4)]);


