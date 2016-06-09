function [dipol] = dipol(Para, Prep, Data)

dipol = cell(1,4);

[grad_H_kx , grad_H_ky] = grad_TB_Liu_TNN_fun(Data.k(1:2,:,1),Para.TB);
mapping = reshape(1:9,[3,3]);

for nll = 1:Para.nr.dipol
        
    hh = Para.dipol_trans(nll,1);
    ee = Para.dipol_trans(nll,2);
    
    if hh <= 3
        d = 0;
    else
        d = 3;
    end
    
    dipol_k = zeros(2, Para.nr.k);
    
    for a = 1:3
        
        for b = 1:3
            
            dipol_k(1,:) = dipol_k(1,:) + grad_H_kx(mapping(a,b),:) .* Prep.CV(:,1,a+d,b+d,hh,ee).';
            dipol_k(2,:) = dipol_k(2,:) + grad_H_ky(mapping(a,b),:) .* Prep.CV(:,1,a+d,b+d,hh,ee).';
                        
        end
        
    end
    
    dipol{1,nll} = Para.vorf.dipol / 1i  * dipol_k ./ repmat( Data.Ek(hh,:,1) - Data.Ek(ee,:,1) ,2,1);
    
end

