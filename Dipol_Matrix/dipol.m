function [dipol] = dipol(Para, Prep, Data)

dipol = cell(6);

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
            
            strAR = ['dA_H_R_',num2str(a),num2str(b),'_rvec.txt'];
            strAI = ['dA_H_I_',num2str(a),num2str(b),'_rvec.txt'];
            
            strBR = ['dB_H_R_',num2str(a),num2str(b),'_rvec.txt'];
            strBI = ['dB_H_I_',num2str(a),num2str(b),'_rvec.txt'];
            
            AR  = load(strAR);
            AI = load(strAI);
            
            BR  = load(strBR);
            BI = load(strBI);
            
            fprintf('a = %d und b = %d',[a,b])
            
            max(abs(real(grad_H_kx(mapping(a,b),:)).' - AR ))
            max(abs(imag(grad_H_kx(mapping(a,b),:)).' - AI ))
            
            max(abs(real(grad_H_ky(mapping(a,b),:)).' - BR ))
            max(abs(imag(grad_H_ky(mapping(a,b),:)).' - BI ))
            
            
        end
        
    end
    
    dipol{n,m} = Para.vorf.dipol / 1i  * dipol_k ./ repmat( Data.Ek(m,:,1) - Data.Ek(n,:,1) ,2,1);
    
end

