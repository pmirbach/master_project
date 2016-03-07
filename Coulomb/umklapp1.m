function [kneu] = umklapp1(Parameter,k,k0)

b1 = Parameter.rezGV(:,1);          % Reziproke Gittervektoren
b2 = Parameter.rezGV(:,2);

b_0 = [0; 0];                       % Same red BZ
b_nn = [1 0 1 -1 0 -1 ; ...         % Next neighbours
    0 1 1 0 -1 -1];
b_snn = [2 2 2 1 0 -1 -2 -2 -2 -1 0 1 ; ... % Second next neighbour
    0 1 2 2 2 1 0 -1 -2 -2 -2 -1];
c = [b_0, b_nn, b_snn];

% All displacement vectors
b = zeros(3, size(c,2));            
for ii = 1:size(c,2)
    b(1:2,ii) = c(1,ii) * b1 + c(2,ii) * b2;
end
 
newks = zeros(3,size(k,2),size(k,3),size(b,2));
newks_diff = zeros(3,size(k,2),size(k,3),size(b,2));
k0m = [k0(1) * ones(1,size(k,2),size(k,3)) ; ...
    k0(2) * ones(1,size(k,2),size(k,3)) ; ...
    zeros(1,size(k,2),size(k,3))];
for ii = 1:size(b,2)        % Add every neighbor vector
    newks(:,:,:,ii) = k + ...
        repmat( b(:,ii) * ones(1,size(k,2)) ,1,1,size(k,3));
    newks_diff(:,:,:,ii) = k - k0m + ...
        repmat( b(:,ii) * ones(1,size(k,2)) ,1,1,size(k,3));
end

difference = sqrt(newks_diff(1,:,:,:).^2 + newks_diff(2,:,:,:).^2); % Norm
% of every new difference vector
min_diff = squeeze(min(difference,[],4));    % Closest k' to k0


kneu = zeros(size(k));
for ii = 1:size(k,2)
    for jj = 1:size(k,3)
        ind = find(difference(1,ii,jj,:) == min_diff(ii,jj));
        if numel(ind) > 1
            ind = ind(1);
        end
        kneu(:,ii,jj) = newks(:,ii,jj,ind);
    end
end
