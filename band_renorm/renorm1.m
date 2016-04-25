function renorm1(Parameter, Ek, V_f, V_h, fk, wk)


[A,B] = meshgrid(1:3,1:3);
c=cat(2,A,B);
ll=reshape(c,[],2);
ll = [ll; ll+3];


Ek_h = Ek;
Ek_f = Ek;
Ek_hf = Ek;

% Parameter.area_BZ

ak = ( fk(1,:) .* wk );
dd = ak * V_f(:,:,1)';
ee = [1 -1 -1];

for ii = 1:3
    
    Ek_f(1,:,1) = Ek_f(1,:,1) + 1 / ( 2 * pi )^2 * ee(ii) * ( fk(ii,:) .* wk ) * V_f(:,:,ii)';
    
end

1


for nk = 1:6
    
    1
    
end