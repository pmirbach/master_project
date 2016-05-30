function [dipol] = dipol(Para, Prep, Data)

dipol = cell(1,4);

transitions = Para.dipol_trans;

[grad_H_kx , grad_H_ky] = grad_TB_Liu_TNN_fun(Data.k(1:2,:,1),Para.TB);
mapping = reshape(1:9,[3,3]);

for nn = 1:size(transitions,1)
    
    m = transitions(nn,1);
    n = transitions(nn,2);
    
    if m <= 3
        d = 0;
    else
        d = 3;
    end
    
    dipol_k = zeros(2, Para.nr.k);
    
    for a = 1:3
        
        for b = 1:3
            
            dipol_k(1,:) = dipol_k(1,:) + grad_H_kx(mapping(a,b),:) .* Prep.CV(:,1,a+d,b+d,m,n).';
            dipol_k(2,:) = dipol_k(2,:) + grad_H_ky(mapping(a,b),:) .* Prep.CV(:,1,a+d,b+d,m,n).';
                        
        end
        
    end
    
    dipol{1,nn} = Para.vorf.dipol / 1i  * dipol_k ./ repmat( Data.Ek(m,:,1) - Data.Ek(n,:,1) ,2,1);
    
end

