function [Ek_hf, Ek_h, Ek_f] = renorm2(Parameter, Ek, V_f, V_h, fk, wk, ll)


Ek_h = Ek(:,:,1);
Ek_f = Ek(:,:,1);
Ek_hf = Ek(:,:,1);


% vorz = [-1 1 1; 1 -1 -1; 1 -1 -1];
% vorz = repmat(vorz,[2,2]);
% 
% for nll = 3 %size(ll,1)
%     
%     l1 = ll(nll,1);
%     l2 = ll(nll,2);
%       
%     for nk = 1:Parameter.nrk
%         
%         V = 0;
%         
%         for nks = 1:Parameter.nrk
%             
%             V = V + fk(l2,nks) * wk(nks) * (-1) * V_f(nk,nks,l2);
%             
%         end
%         
%         Ek_f(l1,nk) = Ek_f(l1,nk) + vorz(l1,l2) * 1 / ( 2 * pi )^2 *  V;
%         
%     end
%     
% end

vorz = [-1 1 1; 1 -1 -1; 1 -1 -1];
vorz = blkdiag(vorz,vorz);

for nll = 1:size(ll,1)
    
    l1 = ll(nll,1);
    l2 = ll(nll,2);
    l3 = ll(nll,3);
   
    Ek_f(l1,:) = Ek_f(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ( fk(l2,:) .* wk / 6 ) * (-1) * V_f(:,:,l3).';
%     Ek_h(l1,:) = Ek_h(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ( fk(l2,:) .* wk / 6 ) * V_h(:,:,l2)';
%     Ek_hf(l1,:) = Ek_hf(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ( fk(l2,:) .* wk / 6 ) * ( V_h(:,:,l2) - V_f(:,:,l2))';
    
end

