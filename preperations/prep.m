function [ CV, CVnoSOC, minq, V_ab_interpl ] = prep( Ctrl, Para, Data, EvnoSOC , Coul_ME)

% Eks = call_eks(Data.Ek);

CV = call_CV( Data.Ev );
CVnoSOC = call_CV( EvnoSOC );

minq = call_minq(Para,Data.k);

V_ab_interpl = get_background_screening( Ctrl , Para , Coul_ME , minq );