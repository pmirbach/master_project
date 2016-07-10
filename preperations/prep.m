function [ CV, CVnoSOC, minq] = prep(Para, Data, EvnoSOC , Coul_ME)

% Eks = call_eks(Data.Ek);

CV = call_CV( Data.Ev );
CVnoSOC = call_CV( EvnoSOC );

minq = call_minq(Para,Data.k);


eps_2 = 0;
eps_3 = 0;
get_background_screening( Para , Coul_ME , eps_2 , eps_3 , max(minq(:)) )