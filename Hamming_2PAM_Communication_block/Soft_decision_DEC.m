%% Haming code soft decision decoding

function decoded_ouptut = Soft_decision_DEC(received_signal)
    code = [0	0	0	0	0	0	0;
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
                1	1	1	1	1	1	1;
end