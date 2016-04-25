function [Ek_hf,Ek_h,Ek_f] = coulomb_3(constAg,Parameter,Data,Prep)

% k = Data.k;

numk = size(Data.k,2);

vorf = constAg.ec^2 / ( 2 * constAg.eps_0);

Ek_h = Data.Ek;
Ek_f = Data.Ek;
Ek_hf = Data.Ek;

renorm_sign = ones(6);
renorm_sign([1,4],[2,3,5,6]) = -1;
renorm_sign([2,3,5,6],[1,4]) = -1;


CV = Prep.CV;
coul_h = Prep.coul_h;


% V_h_test = zeros(1,numk,6);
% V_f_test = zeros(1,numk,6);

ll = [1 1; 1 2];
numll = size(ll,2);

tic

for nk = 1:numk
    
    for nks = 1:numk
        
        for nl = 1:numll
                             
                V_h = 0;
                
                V_f = 0;
                
                for tri = 1:6
                    
                    for a = 1:3
                        
                        for b = 1:3
                            
                            V_h = V_h + ( real( CV(a,a,l1,l1,nk,1) * CV(b,b,l2,l2,nks,tri) ) + ...
                                real( CV(a+3,a+3,l1,l1,nk,1) * CV(b+3,b+3,l2,l2,nks,tri) ) ) * ...
                                coul_h(a,b);
                            
                            V_f = V_f + ( real( CV(a,b,l1,l1,nk,1) * CV(b,a,l2,l2,nks,tri) ) + ...
                                real( CV(a+3,b+3,l1,l1,nk,1) * CV(b+3,a+3,l2,l2,nks,tri) ) ) ...
                                * Prep.coul_intrp{a,b}(Prep.minq(nk,nks,tri));
                            
                        end
                        
                    end
                                        
                end
                
                
                %                                 Ek_h(l1,nk) = Ek_h(l1,nk) + renorm_sign(l1,l2) * ...
                %                                     Data.k(3,nks,1) * vorf * V_h * Data.fk(l2,nks);
                %
                %                                 Ek_f(l1,nk) = Ek_f(l1,nk) + renorm_sign(l1,l2) * ...
                %                                     Data.k(3,nks,1) * vorf * ( - V_f ) * Data.fk(l2,nks);
                %
                %                                 Ek_hf(l1,nk) = Ek_hf(l1,nk) + renorm_sign(l1,l2) * ...
                %                                     Data.k(3,nks,1) * vorf * ( V_h - V_f ) * Data.fk(l2,nks);
                
                
            end
            
            
            %             figure; hold on
            %
            %             for ii = 1:6
            %                 scatter3(Data.k(1,:,ii),Data.k(2,:,ii),V_f_test(1,:,ii)')
            %             end
            
            %             1
        
    end
    
end
toc

