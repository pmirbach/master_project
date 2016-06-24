function [Eks] = call_eks(Ek)

% Übergang ins Elektron-Loch Bild:
% Die Energiebänder werden so verschoben, dass der niedrigste Punkt aller
% Leitungsbänder bei 0 liegt, sowie der oberste Punkte aller Valenzbänder.
% Anschließend wechseln die Valenzbänder das Vorzeichen, was einer
% Berechnung der Löcher entspricht.

band_ind{1} = [1,4];        % Indizes der Valenzbänder
band_ind{2} = [2,3];    % Indizes der Leitungsbänder mit Spin up
band_ind{3} = [5,6];    % Indizes der Leitungsbänder mit Spin down

Eks(band_ind{1},:) = - Ek(band_ind{1},:,1);
Eks(band_ind{2},:) = Ek(band_ind{2},:,1) - min(min(Ek(band_ind{2},:,1)));
Eks(band_ind{3},:) = Ek(band_ind{3},:,1) - min(min(Ek(band_ind{3},:,1)));


% Eks = zeros( size(Ek,1) , size(Ek,2) );
% for ii = band_ind{1}
%     Eks(ii,:) = - Ek(ii,:,1);
% end
% for ii = band_ind{2}
%     Eks(ii,:) = Ek(ii,:,1) - min(min(Ek(band_ind{2},:,1)));
% end

% Tests
if any( Eks( Eks < 0 ) )
    warning('Energybands in electron-hole picture smaller than zero!')
end
lowestE = min(Eks([1,2,4,5],:).');
if any( lowestE( lowestE > 0 ) )
    warning('Valence or lower conduction bands in electron-hole picture do not have min(E) = 0!')
end

