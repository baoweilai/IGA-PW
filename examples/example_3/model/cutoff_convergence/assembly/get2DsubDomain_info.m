function dg_2D_domains_idx = get2DsubDomain_info(mp)
% Build the 2-D subdomain connectivity data.
nx = mp.nx;
ny = mp.ny;
%% For elements from left to right, layer to layer format
dg_2D_domains_idx_left2right = zeros(ny,nx); % The i-th row contains the index along the x-axis

for j = 1:ny
    local_idx = 1;
    for i=1:nx
        e = i + (j-1)*nx;
        dg_2D_domains_idx_left2right(j,local_idx) = e;
        local_idx = local_idx + 1;
    end
end
%% For elements from bottom to top, layer to layer format

dg_2D_domains_idx_bottom2top = zeros(nx,ny); % The i-th row contains the index along the x-axis

for i = 1:nx
    local_idx = 1;
    for j = 1:ny
        e = i + (j-1)*nx;
        dg_2D_domains_idx_bottom2top(i,local_idx) = e;
        local_idx = local_idx + 1;
    end
end

dg_2D_domains_idx.dg_2D_domains_idx_left2right = dg_2D_domains_idx_left2right;
dg_2D_domains_idx.dg_2D_domains_idx_bottom2top = dg_2D_domains_idx_bottom2top;


end
