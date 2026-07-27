function  bnd_dofs_idx = get2DsubDomain_bnd_dof(nurbs_refine_domains,mp)
% Find boundary degrees of freedom for each subdomain.
nx = mp.nx;
ny = mp.ny;
%% For elements on left edge

left_elements_idx = zeros(ny,1);

i = 1;
for j = 1:ny
    e = i + (j-1)*nx;
    left_elements_idx(j) = e;
end

%% For elements on right edge

right_elements_idx = zeros(ny,1);

i = nx;
for j = 1:ny
    e = i + (j-1)*nx;
    right_elements_idx(j) = e;
end

%% For elements on bottom edge

bottom_elements_idx = zeros(nx,1);

j = 1;
for i = 1:nx
    e = i + (j-1)*nx;
    bottom_elements_idx(i) = e;
end

%% For elements on right edge

top_elements_idx = zeros(nx,1);

j = ny;
for i = 1:nx
    e = i + (j-1)*nx;
    top_elements_idx(i) = e;
end
%%  Find the dof indices along the bottom edge
n_dofs_bottom_edge = 0;


for i = 1:nx
    n_dofs_bottom_edge = n_dofs_bottom_edge + nurbs_refine_domains{bottom_elements_idx(i)}.n_dofs_bottom;
end

dofs_idx_bottom = zeros(1,n_dofs_bottom_edge);

local_idx = 1;
for i = 1:nx
    bottom_dofs = nurbs_refine_domains{bottom_elements_idx(i)}.bottom_dofs;
    for d = 1:nurbs_refine_domains{bottom_elements_idx(i)}.n_dofs_bottom
        dofs_idx_bottom(local_idx) = bottom_dofs(d);
        local_idx = local_idx + 1;
    end
end


%%  Find the dof indices along the top edge
n_dofs_top_edge = 0;


for i = 1:nx
    n_dofs_top_edge = n_dofs_top_edge + nurbs_refine_domains{top_elements_idx(i)}.n_dofs_top;
end

dofs_idx_top = zeros(1,n_dofs_top_edge);

local_idx = 1;
for i = 1:nx
    top_dofs = nurbs_refine_domains{top_elements_idx(i)}.top_dofs;
    for d = 1:nurbs_refine_domains{top_elements_idx(i)}.n_dofs_top
        dofs_idx_top(local_idx) = top_dofs(d);
        local_idx = local_idx + 1;
    end
end

%%  Find the dof indices along the left edge
n_dofs_left_edge = 0;


for j = 1:ny
    n_dofs_left_edge = n_dofs_left_edge + nurbs_refine_domains{left_elements_idx(j)}.n_dofs_left;
end

dofs_idx_left = zeros(1,n_dofs_left_edge);

local_idx = 1;
for j = 1:ny
    left_dofs = nurbs_refine_domains{left_elements_idx(j)}.left_dofs;
    for d = 1:nurbs_refine_domains{left_elements_idx(j)}.n_dofs_left
        dofs_idx_left(local_idx) = left_dofs(d);
        local_idx = local_idx + 1;
    end
end

%%  Find the dof indices along the right edge
n_dofs_right_edge = 0;


for j = 1:ny
    n_dofs_right_edge = n_dofs_right_edge + nurbs_refine_domains{right_elements_idx(j)}.n_dofs_right;
end

dofs_idx_right = zeros(1,n_dofs_right_edge);

local_idx = 1;

for j = 1:ny
    right_dofs = nurbs_refine_domains{right_elements_idx(j)}.right_dofs;
    for d = 1:nurbs_refine_domains{right_elements_idx(j)}.n_dofs_right
        dofs_idx_right(local_idx) = right_dofs(d);
        local_idx = local_idx + 1;
    end
end
bnd_dofs_idx = [dofs_idx_bottom,dofs_idx_top,dofs_idx_left,dofs_idx_right];

bnd_dofs_idx = unique(bnd_dofs_idx);


end
