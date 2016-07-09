function [ ev , ef ] = solve_sort_eig( H )

[ef, ev] = eig( H );

[ ev , I ] = sort(diag(real(ev)));
ef = ef(:, I);
