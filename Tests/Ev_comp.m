function [compa] = Ev_comp(Ev1, Ev2)

% Ev: 4-dim
% alpha, lambda, k, tri

s_Ev1 = size(Ev1);
s_Ev2 = size(Ev2);

if ~all( s_Ev1 == s_Ev2 )
    error('Eigenvalues have different sizes!')
end

compa = zeros(s_Ev1(2), s_Ev1(3), s_Ev1(4));
    
for nl = 1:s_Ev1(2)
    for nk = 1:s_Ev1(3)
        for tri = 1:s_Ev1(4)
            
%             diff = Ev1(:,nl,nk,tri) - Ev2(:,nl,nk,tri);
            skal = dot( Ev1(:,nl,nk,tri), Ev2(:,nl,nk,tri) );
            
            phi = rad2deg( atan2( imag(skal), real(skal) ) );
            
            compa(nl,nk,tri) = phi;
            
%             if any(diff)
%                 
%                 if any(diff > 1e-8)
%                     compa(nl,nk,tri) = 1;
%                 end
%                 
%             end
            
        end
    end
end