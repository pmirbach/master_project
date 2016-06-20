function [newvalue, dim_string] = convert_SI_natural(value, dim1, dim2, dim3, dim4, direction)
% [newvalue, dim_string] = convert_SI_natural(value, dim1, dim2, dim3,dim4,
% direction)
% Basis Wechsel Matrizen von SI zu natürlichen Einheiten und zurück
% 3 Inputs ab dem zweiten: Potenzen der Dimensionen der umzurechnenden
%
% Größe: 
%       m^(dim1)    s^(dim2)     kg^(dim3)   A^(dim4) 
%    bzw.
%       c^(dim1)    hquer^(dim2) ev^(dim3)   (4 pi e0)^(dim4)
%
% letzter Input: 
%    string 'SI' für von SI zu natural
%           'nat' für von natural zu SI
%
% siehe http://www.uni-bonn.de/~gansmann/Seite/Uni/Einheiten.pdf
% damit kommt ilans tabelle richtig raus :)

% Basis-Vektoren der natürlichen Einheiten geschrieben in SI-Einheiten:

showOutput = 0;

value_log = log10(value);
c=2.99792458e8; %m/s
hquer = 1.054571628e-34; %J/s = m^2 kg /s
eV = 1.602176487e-19; %J = m^2 kg / s^2
e0 = 8.854187817e-12; %A^2 s^4 kg^-1 m^-2

% Spalten der Matrix: nat->SI: log10(Basis), Potenz der Länge, Potenz der
% Zeit, Potenz der Masse

zehn_vec = [1; 0; 0; 0; 0];
c_vec = [log10(c); 1; -1; 0; 0];
eV_vec = [log10(eV); 2; -2; 1; 0];
hquer_vec = [log10(hquer); 2; -1; 1; 0];
e0_vec = [log10(4*pi*e0); -3; 4; -1; 2];

mat_nat_to_SI = [zehn_vec, c_vec, hquer_vec, eV_vec, e0_vec ];

if strcmp('SI',direction)
   new = mat_nat_to_SI\[value_log,dim1,dim2,dim3,dim4]';
elseif strcmp('nat',direction)
   new = mat_nat_to_SI*[value_log,dim1,dim2,dim3,dim4]';  
else (strcmp(direction,'SI') == 0) & (strcmp(direction,'nat'))
   error('Error: Please type "SI" or "nat" for direction(last Input)!')
end

newvalue = 10^(new(1));
if strcmp('nat',direction)
     dim_string_old = [num2str(value,'%10.3e\n'),' c^',num2str(dim1),' hquer^',num2str(dim2),' eV^',num2str(dim3),' (4pie0)^',num2str(dim4)];
   dim_string = [num2str(newvalue,'%10.3e\n'),' m^',num2str(new(2)),' s^',num2str(new(3)),' kg^',num2str(new(4)),' A^',num2str(new(5))];
elseif strcmp('SI',direction)
   dim_string_old = [num2str(value,'%10.3e\n'),' m^',num2str(dim1),' s^',num2str(dim2),' kg^',num2str(dim3),' A^',num2str(dim4)];
   dim_string = [num2str(newvalue,'%10.3e\n'),' c^',num2str(new(2)),' hquer^',num2str(new(3)),' eV^',num2str(new(4)),' (4pie0)^',num2str(new(5))];
end
if showOutput == 1
    disp(['converted ',10,dim_string_old,'   =>     ',dim_string])
end
