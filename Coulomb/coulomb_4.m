function coulomb_4(constAg,Parameter,Data,Prep)

numk = size(Data.k,2);

vorf = constAg.ec^2 / ( 2 * constAg.eps_0);


[A,B] = meshgrid(1:3,1:3);
c=cat(2,A',B');
ll=reshape(c,[],2);


V_f = zeros(numk,numk,6);
V_h = V_f;

tic

for nll = 1 %:size(ll,1)
    
    l1 = ll(nll,1);
    l2 = ll(nll,2);
    
    for nk = 1:numk
        
        %             V_f_test = zeros(numk,6);
        %             V_h_test = V_f_test;
        %
        %             V_f_o_test = V_f_test;
        %             V_h_o_test = V_f_test;
        
        for nks = 1:numk
            
            for tri = 1:6
                
                for a = 1:3
                    
                    for b = 1:3
                        
                        V_f(nk,nks,tri) = V_f(nk,nks,tri) + ...
                            real( vorf * Prep.CV(a,b,l1,l1,nk,1) * Prep.CV(b,a,l2,l2,nks,tri) * Prep.coul_intrp{a,b}(Prep.minq(nk,nks,tri)) );
                        
                        V_h(nk,nks,tri) = V_h(nk,nks,tri) + ...
                            real( vorf * Prep.CV(a,a,l1,l1,nk,1) * Prep.CV(b,b,l2,l2,nks,tri) * Prep.V_orbital_h(a,b) );
                        
                        
                        
                        %                             V_f_test(nks,tri) = V_f_test(nks,tri) + ...
                        %                                 real( vorf * Prep.CV(a,b,l1,l1,nk,1) * Prep.CV(b,a,l2,l2,nks,tri) * Prep.coul_intrp{a,b}(Prep.minq(nk,nks,tri)) );
                        %
                        %                             V_h_test(nks,tri) = V_h_test(nks,tri) + ...
                        %                                 real( vorf * Prep.CV(a,a,l1,l1,nk,1) * Prep.CV(b,b,l2,l2,nks,tri) * Prep.V_orbital_h(a,b) );
                        %
                        %                             V_f_o_test(nks,tri) = V_f_o_test(nks,tri) + Prep.coul_intrp{a,b}(Prep.minq(nk,nks,tri));
                        %                             V_h_o_test(nks,tri) = V_h_o_test(nks,tri) + Prep.V_orbital_h(a,b);
                        
                        
                        
                    end
                    
                end
                
            end
            
        end
        
        %             figure
        %             hold on
        %             for ii = 1:6
        %                 scatter3(Data.k(1,:,ii),Data.k(2,:,ii),V_f_test(:,ii)')
        %             end
        %
        %             1
        
    end
    
    
end

toc

figure
hold on
for ii = 1:6
    scatter3(Data.k(1,:,ii),Data.k(2,:,ii),V_f(1,:,ii)')
end

% 1
