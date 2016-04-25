function [Ek_hf, Ek_h, Ek_f] = renorm2(Parameter, Ek, V_f, V_h, fk, wk)


[A,B] = meshgrid(1:3,1:3);
c=cat(2,A,B);
ll=reshape(c,[],2);
ll = [ll; ll+3];


Ek_h = Ek(:,:,1);
Ek_f = Ek(:,:,1);
Ek_hf = Ek(:,:,1);

% Parameter.area_BZ

ak = ( fk(1,:) .* wk );
dd = ak * V_f(:,:,1)';
ee = [1 -1 -1];

% for ii = 1:3
%     
%     Ek_f(1,:) = Ek_f(1,:) + 1 / ( 2 * pi )^2 * ee(ii) * ( fk(ii,:) .* wk ) * V_f(:,:,ii)';
%     
% end

vorz = [1 -1 -1 -1 1 1 -1 1 1];
vorz = [vorz, vorz];


for nll = 1:size(ll,1)
    
    l1 = ll(nll,1);
    l2 = ll(nll,2);
   
    Ek_f(l1,:) = Ek_f(l1,:) + 1 / ( 2 * pi )^2 * vorz(nll) * ( fk(l2,:) .* wk ) * (-1) * V_f(:,:,l2)';
    Ek_h(l1,:) = Ek_h(l1,:) + 1 / ( 2 * pi )^2 * vorz(nll) * ( fk(l2,:) .* wk ) * V_h(:,:,l2)';
    Ek_hf(l1,:) = Ek_hf(l1,:) + 1 / ( 2 * pi )^2 * vorz(nll) * ( fk(l2,:) .* wk ) * ( V_h(:,:,l2) - V_f(:,:,l2))';
    
end


% for nk = 1:6
%     
%     1
%     
% end