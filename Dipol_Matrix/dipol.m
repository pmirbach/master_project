function [dipol] = dipol(Parameter, Data)


dipol = cell(6);

transitions = [1, 2 ; 1 , 3 ; 2, 1 ; 3 , 1 ];
transitions = [transitions ; transitions + 3];

for nn = 1:size(transitions,1)
    
    m = transitions(nn,1);
    n = transitions(nn,2);
    
    dipol_k = zeros(2, size(Data.k,1));
    
    for ii = 1:size(Data.k,2)
        prefix = 1.602e4 / ( 1i * ( Data.Ek(m,ii) - Data.Ek(n,ii) ) );
        [grad_H_kx , grad_H_ky] = ...
            grad_TB_Liu_TNN_fun(Data.k(1:2,ii,1) , Parameter);
%         B =  transpose(Data.Ev(:,n,ii) * Data.Ev(:,m,ii)');
        transform_tensor = conj(Data.Ev(:,n,ii))*Data.Ev(:,m,ii).';
               
        summand_kx = blkdiag(grad_H_kx,grad_H_kx) .* transform_tensor;        
        sum_kx = sum(summand_kx(:));     
        dipol_k(1,ii) = prefix * sum_kx ;
           
        summand_ky = blkdiag(grad_H_ky,grad_H_ky) .* transform_tensor;       
        sum_ky = sum(summand_ky(:));
        dipol_k(2,ii) = prefix * sum_ky ;
    end
    
    dipol{n,m} = dipol_k;
    
%     dipol_plot(nn,:) = abs(1 / sqrt(2) * ...
%         (dipol{n,m}(1,:) - 1i * dipol{n,m}(2,:)));
    
end
