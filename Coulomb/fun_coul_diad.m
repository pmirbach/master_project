function [coul_diad_h, coul_diad_f] = fun_coul_diad(Ev_k, Ev_ks, combi)

% coul_diad_h = zeros(3,3,6,2);
% coul_diad_f = zeros(3,3,6,2);
%
% Evk(:,:,1) = Ev_k(1:3,1:3);
% Evk(:,:,2) = Ev_k(4:6,4:6);
%
% Evks(:,:,:,1) = squeeze( Ev_ks(1:3,1:3,:,:) );
% Evks(:,:,:,2) = squeeze( Ev_ks(4:6,4:6,:,:) );

% for a = 1:3
%
%     for b = 1:3
%
%         for ni = 1:6
%
%             for s = 1:2
%
%
%                 coul_diad_h(a,b,ni,s) = conj(Evk(a,combi(1),s)) * ...
%                     conj(Evks(b,combi(2),ni,s)) * ...
%                     Evks(b,combi(3),ni,s) * Evk(a,combi(4),s);
%
%                 coul_diad_f(a,b,ni,s) = conj(Evk(a,combi(1),s)) *...
%                     conj(Evks(b,combi(2),ni,s)) * ...
%                     Evk(b,combi(3),s) * Evks(a,combi(4),ni,s);
%
%             end
%
%         end
%
%     end
%
% end

coul_diad_h = zeros(6,6,6);
coul_diad_f = zeros(6,6,6);

for a = 1:6
    
    for b = 1:6
        
        for ni = 1:6
            
            coul_diad_h(a,b,ni) = conj(Ev_k(a,combi(1))) * ...
                Ev_k(a,combi(4)) * conj(Ev_ks(b,combi(2),ni)) * ...
                Ev_ks(b,combi(3),ni) ;
            
            coul_diad_f(a,b,ni) = conj(Ev_k(a,combi(1))) *...
                conj(Ev_ks(b,combi(2),ni)) * ...
                Ev_k(b,combi(3)) * Ev_ks(a,combi(4),ni) ;
        end
        
    end
    
end
