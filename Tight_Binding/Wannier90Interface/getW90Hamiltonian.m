function [ W90Hamiltonian , varargout ]= getW90Hamiltonian(W90Data, kk)

% reconstruct H(k)_mn for every height
%
%   H(k)_mn = SUM_R_i 1/degeneracy(R_i) H(R_i)_mn * exp(1i * k * R_i)
%
%       with: R_i = hamiltonianRaw(i, 1:3) * rlat
%


W90HamiltonianTmp = zeros(size(kk, 1), W90Data(1).numWann, W90Data(1).numWann);
W90HamiltonianGradKx = zeros(size(kk, 1), W90Data(1).numWann, W90Data(1).numWann);
W90HamiltonianGradKy = zeros(size(kk, 1), W90Data(1).numWann, W90Data(1).numWann);

indDegen = reshape( repmat( 1:W90Data.numWSGridPoints , 9 , 1) , 1 , [] );

for nn = 1 : length(W90Data.hamiltonianRaw)
    
    Hij = ( W90Data.hamiltonianRaw(nn, 6)  ...
        * exp(1i * 2*pi * kk(:, 1:2) * W90Data.hamiltonianRaw(nn, 1:2)') ...
        * 1 / W90Data.degeneracy( indDegen(nn) ) );
    
    W90HamiltonianTmp(:, W90Data.hamiltonianRaw(nn, 4), W90Data.hamiltonianRaw(nn, 5)) = ...
        W90HamiltonianTmp(:, W90Data.hamiltonianRaw(nn, 4), W90Data.hamiltonianRaw(nn, 5)) ...
        + Hij;
    
    if nargout > 1
        
        W90HamiltonianGradKx(:, W90Data.hamiltonianRaw(nn, 4), W90Data.hamiltonianRaw(nn, 5)) = ...
            W90HamiltonianGradKx(:, W90Data.hamiltonianRaw(nn, 4), W90Data.hamiltonianRaw(nn, 5)) ...
            + ( 1i* ( W90Data.hamiltonianRaw(nn, 1:3) * W90Data.rLat(:,1) ) ...   % i*Rx  out  of  dx e(ikR) = iRx * e(ikR)
            * Hij );
        
        W90HamiltonianGradKy(:, W90Data.hamiltonianRaw(nn, 4), W90Data.hamiltonianRaw(nn, 5)) = ...
            W90HamiltonianGradKy(:, W90Data.hamiltonianRaw(nn, 4), W90Data.hamiltonianRaw(nn, 5)) ...
            + ( 1i* ( W90Data.hamiltonianRaw(nn, 1:3) * W90Data.rLat(:,2) ) ...   % i*Ry  out  of  dy e(ikR) = iRy * e(ikR)
            * Hij );
        
    end
    
end

varargout{1} = W90HamiltonianGradKx;
varargout{2} = W90HamiltonianGradKy;



% Include spin (... or not?)
if(~strcmp(W90Data.SOCSettings.type, 'none'))
    
    W90Hamiltonian = zeros(size(kk, 1), W90Data(1).numWann*2, W90Data(1).numWann*2);
    
    for kkii = 1 : size(kk, 1)
        
        % define lambda
        if(strcmp(W90Data.SOCSettings.type, 'first'))
            lambda = W90Data.SOCSettings.Delta / 2;
        end
        if(strcmp(W90Data.SOCSettings.type, 'second'))
            xi = W90Data.SOCSettings.BG - W90Data.SOCSettings.E0p1;
            lambda = sqrt( xi * ( W90Data.SOCSettings.Delta + xi ) ) - xi ;
        end
        
        %  damped, k-dependent SOC
        if(W90Data.SOCSettings.dampedLambda0 == 1)
            ee      = exp(1);
            kkiiLat = W90Data(1).kLat' *  kk(kkii, 1:3)';
            KKLat   = W90Data(1).kLat' * [1/3, 1/3, 0]';
            kkAct   = norm(kkiiLat - KKLat);
            KK      = norm(KKLat);
            lambda  = lambda * ee * (1 - kkAct / KK)^2 * exp( -(1 - kkAct / KK)^2 );
        end
        
        % get k-dependent W90 matrix
        HHTmp(1:W90Data(1).numWann, ...
            1:W90Data(1).numWann) = ...
            W90HamiltonianTmp(kkii, ...
            1:W90Data(1).numWann, ...
            1:W90Data(1).numWann);
        
        % coupling matrices
        %   first order
        SOC1 = [  0,   0,   0;
            0,   0, +1i;
            0, -1i,   0];
        
        SOC1Up   = -1 * SOC1 * lambda;
        SOC1Down =      SOC1 * lambda;
        
        %   spin up
        HHTmpUp   = HHTmp + SOC1Up;
        %   spin down
        HHTmpDown = HHTmp + SOC1Down;
        
        %   second order
        if(strcmp(W90Data.SOCSettings.type, 'second'))
            
            a = 3/2 * lambda^2 * 1 / W90Data.SOCSettings.E0p1;
            b =       lambda^2 * 1 / W90Data.SOCSettings.Em2m1;
            c =       lambda^2 * 1 / (W90Data.SOCSettings.E0p1 - W90Data.SOCSettings.BG);
            d = 3/2 * lambda^2 * 1 / W90Data.SOCSettings.E0m1;
            
            SOC2Up = [  2*a, 0,     0;
                0,   b,    1i*b;
                0,  -1i*b,  b     ];
            
            SOC2Down = [  2*d, 0,     0;
                0,   c,    -1i*c;
                0,   1i*c, c     ];
            
            HHTmpUp   = HHTmpUp   + 1/2 * SOC2Up;
            HHTmpDown = HHTmpDown + 1/2 * SOC2Down;
            
        end
        
        % complete
        W90Hamiltonian(kkii, ...
            1:2*W90Data(1).numWann, ...
            1:2*W90Data(1).numWann) = ...
            [HHTmpUp,  zeros(W90Data(1).numWann);
            zeros(W90Data(1).numWann), HHTmpDown];
        
    end
    
    % No spin ...
else
    W90Hamiltonian = W90HamiltonianTmp;
end

end
