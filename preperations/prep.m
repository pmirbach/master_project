function [Eks, CV, CV2, minq, V_orbital_h] = prep(Parameter, Data)


Eks = call_eks(Data.Ek);

CV = calc_CV( Data.Ev );
CV2 = calc_CV2( Data.Ev );

% minq = call_minq(Parameter,Data.k);

minq = call_minq3(Parameter,Data.k);

% [coul_intrp] = call_coul_intrp(constAg, Parameter);

[V_orbital_h] = fun_coul_orbital_hartree(Parameter.coul_screened);
