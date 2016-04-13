function [Ek_hf,Ek_h,Ek_f] = coulomb_2(constAg,Parameter,Data,Prep)

% k = Data.k;

numk = size(k,2);

vorf = constAg.ec^2 / ( 2 * constAg.eps_0);

Ek_h = Data.Ek;
Ek_f = Data.Ek;
Ek_hf = Data.Ek;

renorm_sign = ones(6);
renorm_sign([1,4],[2,3,5,6]) = -1;
renorm_sign([2,3,5,6],[1,4]) = -1;

% tic

for l1 = 1:6
    
    for l2 = 1:6
        
        for nk = 1:numk
            
            for nks = 1:numk
                
                for tri = 1:6
                    
                    do
                    
                end
                
            end
            
        end
        
    end
    
end
        


for nk = 1:size(Data.k,2)
    
    disp(nk)
    
    % Neue k nach Umklapp Prozess
    k_shift = umklapp1(Parameter, Data.k(1:2,:,:), Data.k(1:2,nk,1));
    
    for nks = 1:size(Data.k,2)
                
        q_v = repmat(Data.k(1:2,nk,1),1,1,6) - k_shift(1:2,nks,:);
        q = squeeze(sqrt( q_v(1,1,:).^2 + q_v(2,1,:).^2 ));
        
%         [V_orbital_f] = fun_coul_orbital_fock(q, ...
%             Parameter.coul_screened, Parameter.coul_kappa );
        [V_orbital_f] = fun_coul_orbital_fock2(q, ...
            Parameter.coul_screened, Parameter.coul_kappa);
                      
        for l1 = 1:6            % Zu renormierendes Band
            
            for l2 = 1:6        % Andere Bänder   

%                 [coul_diad_h] = ...
%                     coul_hartree(Data.Ev(:,:,nk,1), Data.Ev(:,:,nks,:), l1, l2);
%                 
                for ntri = 1:6
                    
                    coul_diad_h(:,:,ntri) = diag(CV(:,:,l1,l1,nk,1)) * ...
                        diag(CV(:,:,l2,l2,nks,ntri))';
                    coul_diad_f(:,:,ntri) = CV(:,:,l1,l1,nk,1) .* ...
                        CV(:,:,l2,l2,nks,ntri).';
                    
                end
                
%                 [coul_diad_f] = ...
%                     coul_fock(Data.Ev(:,:,nk,1), Data.Ev(:,:,nks,:), l1, l2);
                
                               
                V_h = real( sum( sum( sum( coul_diad_h .* V_orbital_h ) ) ) ); 
                V_f = real( sum( sum( sum( coul_diad_f .* V_orbital_f ) ) ) );
                          
                Ek_h(l1,nk) = Ek_h(l1,nk) + renorm_sign(l1,l2) * ...
                    Data.k(3,nks,1) * vorf * V_h * Data.fk(l2,nks);
                Ek_f(l1,nk) = Ek_f(l1,nk) + renorm_sign(l1,l2) * ...
                    Data.k(3,nks,1) * vorf * ( - V_f ) * Data.fk(l2,nks);
                Ek_hf(l1,nk) = Ek_hf(l1,nk) + renorm_sign(l1,l2) * ...
                    Data.k(3,nks,1) * vorf * ( V_h - V_f ) * Data.fk(l2,nks);
                       
            end
            
        end  
                
              
    end

end

% toc




