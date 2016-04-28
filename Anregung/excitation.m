function fk = excitation(Ctrl,constAg,Para,k,Eks)

T = Ctrl.temperature;                       % Temperatur
D_0 = Ctrl.carrier_density * 1e-14;         % Anregungsdichte in 1/nm^2
tol = Ctrl.carrier_density_tol * 1e-14;     % Toleranz

fk = zeros(size(Eks));

% Zunächst müssen für Elektronen und Löcher die chemischen Potentiale mu so
% bestimmt werden, dass ihre Anregungsdichten D0 betragen.

for ii = 1:2
    ind = Para.TB_ind{ii};
    Nr_ind = size(ind,2);
    
    mu = [-1000, 0, 1000];      % Grenzen zum Suchen - Verbessern!!!
    
    diff(2) = D_0;
    
    lauf = 0;
    
    while abs(diff(2)) > tol
        mu(2) = (mu(1) + mu(3)) / 2;
        
        W_k = zeros(size(ind,2),size(Eks,2),3);
        D = zeros(size(ind,2),3);
        
        for jj = 1:Nr_ind
            for kk = 1:3
                W_k(jj,:,kk) = 1 ./ ( exp( ( Eks(ind(jj),:) ...
                    - mu(kk) ) ./ ( constAg.k_B * T ) ) + 1 );
                D(jj,kk) = 1 / (2 * pi)^2 * sum( W_k(jj,:,kk) .* k(3,:));
            end
        end
        diff = sum(D) - D_0;
        
        if diff(3) * diff(2) > 0
            mu(3) = mu(2);
        elseif diff(1) * diff(2) > 0
            mu(1) = mu(2);
        else
            error('doof')
        end
        
        lauf = lauf + 1;
        
    end
        
    fk(ind,:) = W_k(:,:,2);
    
end