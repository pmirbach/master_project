function [dipol] = dipol(Parameter, Prep, Data)


dipol = cell(6);


transitions = Parameter.dipol_trans;

for nn = 1:size(transitions,1)
    
    m = transitions(nn,1);
    n = transitions(nn,2);
    
    dipol_k = zeros(2, Parameter.nrk);
    
    for nk = 1:Parameter.nrk
        prefix = 1.602e4 / ( 1i * ( Data.Ek(m,nk) - Data.Ek(n,nk) ) );
        [grad_H_kx , grad_H_ky] = grad_TB_Liu_TNN_fun(Data.k(1:2,nk,1) , Parameter);
%         B =  transpose(Data.Ev(:,n,ii) * Data.Ev(:,m,ii)');
%         transform_tensor = conj(Data.Ev(:,n,nk))*Data.Ev(:,m,nk).';
                 
        dipol_k(1,nk) = prefix * sum( sum( blkdiag(grad_H_kx,grad_H_kx) .* Prep.CV(:,:,n,m,nk,1) ) );
        dipol_k(2,nk) = prefix * sum( sum( blkdiag(grad_H_ky,grad_H_ky) .* Prep.CV(:,:,n,m,nk,1) ) );

    end
    
    dipol{n,m} = dipol_k;
    
%     dipol_plot(nn,:) = abs(1 / sqrt(2) * ...
%         (dipol{n,m}(1,:) - 1i * dipol{n,m}(2,:)));
    
end
