function [B,B_integ] = flaecheninhalt(Parameter,k)

b = norm(Parameter.symmpts{2}(:,2));    % Kantenlänge der BZ    
B = b^2 * 3/2 *sqrt(3);                 % Flächeninhalt der BZ


B_integ = 0;
for ii = 1:size(k,2)
    B_integ = B_integ + k(3,ii);    % Gewicht mal kleines 6-Eck
end
