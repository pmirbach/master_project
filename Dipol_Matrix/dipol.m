function [dipol] = dipol( Para , Prep , Data )

dipol = cell(Para.nr.dipol,2);

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
            
            dipol_k(1,:) = dipol_k(1,:) + Prep.H_grad_kx(mapping(a,b),:) .* Prep.CV(:,1,a+d,b+d,hh,ee).';
            dipol_k(2,:) = dipol_k(2,:) + Prep.H_grad_ky(mapping(a,b),:) .* Prep.CV(:,1,a+d,b+d,hh,ee).';
                        
        end
        
    end
    dipol_k = Para.vorf.dipol / 1i  * dipol_k ./ repmat( Data.Ek(hh,:,1) + Data.Ek(ee,:,1) + Data.EGap ,2,1);
    
    dipol{nll,1} = 1 / sqrt(2) * abs( dipol_k(1,:) + 1i * dipol_k(2,:) ).'; 
    dipol{nll,2} = 1 / sqrt(2) * abs( dipol_k(1,:) - 1i * dipol_k(2,:) ).'; 
     
end

