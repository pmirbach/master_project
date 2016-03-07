function fk = excitation(Ctrl,constAg,k,Ek)

T = Ctrl.temperature;                       % Temperatur
D_0 = Ctrl.carrier_density * 1e-14;         % Anregungsdichte in 1/nm^2
tol = Ctrl.carrier_density_tol * 1e-14;     % Toleranz

fk = zeros(size(Ek));

% Die Energiebänder werden so verschoben, dass der niedrigste Punkt aller
% Leitungsbänder bei 0 liegt, sowie der oberste Punkte aller Valenzbänder.
% Anschließend wechseln die Valenzbänder das Vorzeichen, was einer
% Berechnung der Löcher entspricht.

band_ind{1} = [1,4];        % Indizes der Valenzbänder
band_ind{2} = [2,3,5,6];    % Indizes der Leitungsbänder

Ek_shift(band_ind{1},:) = ...
    -(Ek(band_ind{1},:) - max(max(Ek(band_ind{1},:))));
Ek_shift(band_ind{2},:) = ...
    Ek(band_ind{2},:) - min(min(Ek(band_ind{2},:)));

% Zunächst müssen für Elektronen und Löcher die chemischen Potentiale mu so
% bestimmt werden, dass ihre Anregungsdichten D0 betragen.

for ii = 1:2
    ind = band_ind{ii};
    
    mu = [-1000, 0, 1000];      % Grenzen zum Suchen - Verbessern!!!
    
    diff(2) = D_0;
    
    while abs(diff(2)) > tol
        mu(2) = (mu(1) + mu(3)) / 2;
        
        W_k = zeros(size(ind,2),size(Ek,2),3);
        D = zeros(size(ind,2),3);
        
        for jj = 1:size(ind,2)
            for kk = 1:3
                W_k(jj,:,kk) = 1 ./ ( exp( ( Ek_shift(ind(jj),:) ...
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
        
    end
        
    fk(ind,:) = W_k(:,:,2);
    
end