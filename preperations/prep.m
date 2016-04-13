function [Eks, CV, minq, coul_intrp] = prep(constAg, Parameter, Data)

Eks = call_eks(Data.Ek);

CV = calc_CV( Data.Ev );

minq = call_minq(Parameter,Data.k);

[coul_intrp] = call_coul_intrp(constAg, Parameter);

[V_orbital_h] = fun_coul_orbital_hartree(Parameter.coul_screened);