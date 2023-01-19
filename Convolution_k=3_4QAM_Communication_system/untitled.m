clear
clc
close

t1 = randi([0, 1], [2,100]);

rng('default')
txSig = qammod(t1,4, ...
    InputType='bit', ...
    UnitAveragePower=true);


Esym = sum(abs(txSig).^2)/(length(txSig));
noise_v = 2;
Eb_No = Esym/2/(2*noise_v^2);

rxSig = awgn(txSig, 2*Eb_No, 'measured', 'linear');
rng('default')
my_data=AWGN_Channel(txSig, noise_v);

rxSig - my_data

plot(abs(rxSig), abs(my_data))