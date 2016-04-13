function [Eks] = call_eks(Ek)

% Die Energiebänder werden so verschoben, dass der niedrigste Punkt aller
% Leitungsbänder bei 0 liegt, sowie der oberste Punkte aller Valenzbänder.
% Anschließend wechseln die Valenzbänder das Vorzeichen, was einer
% Berechnung der Löcher entspricht.

band_ind{1} = [1,4];        % Indizes der Valenzbänder
band_ind{2} = [2,3,5,6];    % Indizes der Leitungsbänder

Eks(band_ind{1},:) = - Ek(band_ind{1},:,1);
Eks(band_ind{2},:) = Ek(band_ind{2},:,1) - min(min(Ek(band_ind{2},:,1)));