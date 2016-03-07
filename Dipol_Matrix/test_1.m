clear variables 
close all
clc

B = cell(6);

a = [1 2;3 4; 5 6];

for ii = 1:size(a,1)
    a(ii,:);
    
    B{a(ii,1), a(ii,2)} = 1;
   
end


uebergaenge = [1, 2 ; 1 , 3 ; 2, 3 ; 2, 1 ; 3 , 1 ; 3 , 2];
uebergaenge = [uebergaenge ; uebergaenge + 3];