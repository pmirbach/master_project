function [k, path_tick] = k_path(symm_points, path, gen)

n_path = numel(path);
ind = zeros(1,n_path);
gew = zeros(1,n_path-1);
kx_cell = cell(1,n_path-1);
ky_cell = cell(1,n_path-1);

for ii = 1:n_path
    ind(ii) = find(ismember(symm_points{1},path(ii)));
end

for ii = 1:n_path-1
    gew(ii) = ( symm_points{2}(1,ind(ii)) - symm_points{2}(1,ind(ii+1)) )^2 ...
    + ( symm_points{2}(2,ind(ii)) - symm_points{2}(2,ind(ii+1)) )^2;
end
gew = round(gen * gew / max(gew));

for ii = 1:n_path-1
    kx_cell{ii} = linspace(symm_points{2}(1,ind(ii)),...
        symm_points{2}(1,ind(ii+1)), gew(ii));
    ky_cell{ii} = linspace(symm_points{2}(2,ind(ii)),...
        symm_points{2}(2,ind(ii+1)), gew(ii));
    if ii > 1
        kx_cell{ii}(1) = [];
        ky_cell{ii}(1) = [];
    end
end

k = [cell2mat(kx_cell) ; cell2mat(ky_cell)];

path_tick = zeros(1,n_path);
for ii = 2:n_path
    path_tick(ii) = sum(gew(1:ii-1)) - (ii-1);
end



