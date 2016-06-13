function [V_rabi_fock] = coulomb_rabi_f(Ctrl, Para, Prep, Ev)

V_rabi_fock = zeros( Para.nr.k, Para.nr.k, size(Para.dipol_trans,1) );
V_rabi_fock2 = zeros( Para.nr.k, Para.nr.k, size(Para.dipol_trans,1) );
V_rabi_fock3 = zeros( Para.nr.k, Para.nr.k, size(Para.dipol_trans,1) );



para_map = [1 2 3 ; 2 4 5 ; 3 5 6];

for nll = 1:Para.nr.dipol
    
    hh = Para.dipol_trans(nll,1);
    ee = Para.dipol_trans(nll,2);
    
    if hh <= 3
        d = 0;
    else
        d = 3;
    end
    
    for a = 1:3
        
        for b = 1:3
            
            for tri = 1:6
                
%                 V_rabi_fock(:,:,nll) = V_rabi_fock(:,:,nll) + ...
%                     ( Prep.CV(:,1,a+d,b+d,ee,hh) * Prep.CV(:,tri,b+d,a+d,hh,ee).' ) .* ...
%                     final_coul_scr(Prep.minq(:,:,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                
                V_rabi_fock(:,:,nll) = V_rabi_fock(:,:,nll) + ...
                    ( Prep.CV_noSOC(:,1,a,b,ee-d,hh-d) * Prep.CV_noSOC(:,tri,b,a,hh-d,ee-d).' ) .* ...
                    final_coul_scr(Prep.minq(:,:,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                
                %                 V_rabi_fock(:,:,nll) = V_rabi_fock(:,:,nll) + ...
                %                     final_coul_scr(Prep.minq(:,:,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                
                
                %                 for nk = 1:Para.nr.k
                %
                %                     for nks = 1:Para.nr.k
                % %
                % %                         V_rabi_fock2(nk,nks,nll) = V_rabi_fock2(nk,nks,nll) + ...
                % %                             Prep.CV(nk,1,a+d,b+d,ee,hh) * Prep.CV(nks,tri,b+d,a+d,hh,ee) .* ...
                % %                             final_coul_scr(Prep.minq(nk,nks,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                % %
%                         V_rabi_fock3(nk,nks,nll) = V_rabi_fock3(nk,nks,nll) + ...
%                             conj( Ev(a,ee,nk,1) ) * conj( Ev(b,hh,nks,tri) ) * Ev(b,hh,nk,1) * Ev(a,ee,nks,tri) * ...
%                             final_coul_scr(Prep.minq(nk,nks,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
%                         
%                     end
%                     
%                 end
                
                
            end
            
        end
        
    end
    
end
% % max(max(real(V_rabi_fock(:,:,1)-V_rabi_fock2(:,:,1))))
% % max(max(imag(V_rabi_fock(:,:,1)-V_rabi_fock2(:,:,1))))
% % 
% % max(max(real(V_rabi_fock(:,:,1)-V_rabi_fock3(:,:,1))))
% % max(max(imag(V_rabi_fock(:,:,1)-V_rabi_fock3(:,:,1))))

V_rabi_fock = Para.vorf.coul * abs( V_rabi_fock );
% V_rabi_fock = abs( V_rabi_fock );
