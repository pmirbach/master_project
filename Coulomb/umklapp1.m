function [kneu] = umklapp1(Parameter,k,k0)

b1 = Parameter.rezGV(:,1);          % Reziproke Gittervektoren
b2 = Parameter.rezGV(:,2);

% Vielleicht?
c = [0 1 0 1 -1 0 -1 2 2 2 1 0 -1 -2 -2 -2 -1 0 1;
     0 0 1 1 0 -1 -1 0 1 2 2 2 1 0 -1 -2 -2 -2 -1];
 
% Definition der einstufige reziproken Gittervektoren
b(:,1) = 0 * b1;
b(:,2) = b1;
b(:,3) = b2;
b(:,4) = b1 + b2;
b(:,5) = -b1;
b(:,6) = -b2;
b(:,7) = -(b1 + b2);
% Definition der zweistufige reziproken Gittervektoren
b(:,8) = 2 * b1;
b(:,9) = 2 * b1 + b2;
b(:,10) = 2 * b1 + 2 * b2;
b(:,11) = b1 + 2 * b2;
b(:,12) = 2 * b2;
b(:,13) = -b1 + b2;
b(:,14) = -2 * b1;
b(:,15) = -2 * b1 - b2;
b(:,16) = -2 * b1 - 2 * b2;
b(:,17) = -b1 - 2 * b2;
b(:,18) = -2 * b2;
b(:,19) = b1 - b2;

b = [b ; zeros(1,size(b,2))];

% bm = zeros(3,size(k,2),size(b,2));
newks1 = zeros(3,size(k,2),size(k,3),size(b,2));
newks = zeros(3,size(k,2),size(k,3),size(b,2));
k0m = [k0(1) * ones(1,size(k,2),size(k,3)) ; ...
    k0(2) * ones(1,size(k,2),size(k,3)) ; ...
    zeros(1,size(k,2),size(k,3))];
for ii = 1:size(b,2)
    %     bm(:,:,ii) = repmat( b(:,ii) * ones(1,size(k,2)) ,1,1,size(k,3));
    %     BZ(:,:,ii) = corners + bm(:,:,ii);
    newks1(:,:,:,ii) = k + ...
        repmat( b(:,ii) * ones(1,size(k,2)) ,1,1,size(k,3));
    newks(:,:,:,ii) = k - k0m + ...
        repmat( b(:,ii) * ones(1,size(k,2)) ,1,1,size(k,3));
    %     doof = sum(newks(1:2,:,:,:));
end

doof = sqrt(newks(1,:,:,:).^2 + newks(2,:,:,:).^2);
c = squeeze(min(doof,[],4));


kneu = zeros(size(k));
for ii = 1:size(k,2)
    for jj = 1:size(k,3)
        ind = find(doof(1,ii,jj,:) == c(ii,jj));
        if numel(ind) > 1
            ind = ind(1);
        end
        kneu(:,ii,jj) = newks1(:,ii,jj,ind);
    end
end
