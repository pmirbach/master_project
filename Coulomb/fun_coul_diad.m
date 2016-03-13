function [coul_diad_h, coul_diad_f] = fun_coul_diad(Ev_k, Ev_ks, combi)


coul_diad_h = ( conj(Ev_k(:,combi(1))) .* Ev_k(:,combi(4)) ) * ...
    ( conj(Ev_ks(:,combi(2))) .* Ev_ks(:,combi(3)) )';
% 
% coul_diad_f = ( conj(Ev_k(:,combi(1))) .* Ev_ks(:,combi(4)) ) * ...
%     ( conj(Ev_ks(:,combi(2))) .* Ev_k(:,combi(3)) )';


% coul_diad_h = ( conj(Ev_k(:,combi(1))) * conj(Ev_ks(:,combi(2)))' ) .* ...
%     ( Ev_ks(:,combi(3)) * Ev_k(:,combi(4))' );

coul_diad_f = ( conj(Ev_k(:,combi(1))) * conj(Ev_ks(:,combi(2)))' ) .* ...
    ( Ev_k(:,combi(3)) * Ev_ks(:,combi(4))' );
