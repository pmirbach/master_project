function [Eks, CV, CVnoSOC, minq] = prep(Para, Data, EvnoSOC)

Eks = call_eks(Data.Ek);

CV = call_CV( Data.Ev );
CVnoSOC = call_CV( EvnoSOC );

minq = call_minq(Para,Data.k);