function uh_at_uv   = reference_uh(nurbs_refine,uh,u,v)
% Evaluate the reference solution.

% This script is used to compute the reference solution uh at the
% parametric point (u,v)

pu   =  nurbs_refine.pu;
pv   =  nurbs_refine.pv;

Ubar = nurbs_refine.Ubar;
Vbar = nurbs_refine.Vbar;

m = nurbs_refine.m;

i = findspan(Ubar,pu,u);
j = findspan(Vbar,pv,v);

n_ele_dofs = (pu+1)*(pv+1);

row = zeros(n_ele_dofs,1);
k = 1;

for j1 = (j-pv):j
    for i1 = (i-pu):i
         row(k) = i1 + (j1-1)*m;
         k = k+1;
    end
end

Nu = bsplinebasis(Ubar,pu,u); % A column vector
Nv = bsplinebasis(Vbar,pv,v); % A column vector

N = Nu*(Nv'); N = N(:);

uh_at_uv  = uh(row)'*N;
end
