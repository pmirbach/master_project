function [t_e, x_e] = dgl_rkutta4(t,x,delta,dgl_str)
%%

dgl = str2func(dgl_str);

k1 = dgl(t,x);
k2 = dgl(t + delta/2 , x + delta/2 * k1);
k3 = dgl(t + delta/2 , x + delta/2 * k2);
k4 = dgl(t + delta , x + delta * k3);

x_e = x + delta/6 * (k1 + 2*k2 + 2*k3 + k4);
t_e = t + delta;