function [B_integ] = flaecheninhalt(Para,k)


B_integ = 0;
for ii = 1:size(k,2)
    B_integ = B_integ + k(3,ii);    % Gewicht mal kleines 6-Eck
end
B_integ = Para.BZsmall.area * B_integ;