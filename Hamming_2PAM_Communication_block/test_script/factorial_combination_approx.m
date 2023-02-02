%% x! and (N, r) approximation
clc
clear
close

x =   0:30;
x2 = 0:0.1:30;
x3 = 0:30;
x4 = 0:0.1:30;

poisson = zeros(length(x), 1);
gaussian = zeros(length(x2), 1);
combination = zeros(length(x3), 1);
com_approx = zeros(length(x4), 1);
index(1:4, 1) = 1;


for i = x
    poisson(index(1),1) = (exp(-15)*15.^(i))/(factorial(i));
    index(1) =  index(1) + 1;
end

for j = x2
    gaussian(index(2), 1) = 1/sqrt(2*pi*15)*exp(-((j-15)^2/30));
    index(2) = index(2) + 1;
end

for k = x3
    combination(index(3), 1) = log(factorial(30)/(factorial(30-k)*factorial(k)))
    index(3) = index(3) + 1;
end

for z = x4
    com_approx(index(4), 1) = (30-z)*log(30/(30-z))+z*log(30/z) %- 1/2*log(2*pi*z*(30-z)/30);
    index(4) = index(4) + 1;
end

figure(1),plot(x, poisson, 'ro')
hold on
plot(x2, gaussian, '--')
axis([0 30 0 0.12])
legend('poisson', 'gaussian')

figure(2), plot(x3, combination, 'bx');
hold on
plot(x4, com_approx);
legend('combination', 'approximation')

