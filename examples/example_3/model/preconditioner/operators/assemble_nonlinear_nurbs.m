function N_in = assemble_nonlinear_nurbs( ...
nurbs_original, nurbs_refine, rho_samples, n_gp)

% Assemble the nonlinear NURBS potential matrix.
% Read the refined mesh and sampled density data.
ConPts_o = nurbs_original.ConPts;

Element    = nurbs_refine.Element;
knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
Coordinate = nurbs_refine.Coordinate;
NoEs       = nurbs_refine.NoEs;
n_dofs     = nurbs_refine.n_dofs;
pu         = nurbs_refine.pu;
pv         = nurbs_refine.pv;

rho_samples = rho_samples(:);
if numel(rho_samples) ~= NoEs * n_gp * n_gp
    error('assemble_nonlinear_nurbs: rho_samples length mismatch.');
end

% Prepare quadrature, basis sizes, and sparse matrix triplets.
Fhat = @(x,a,b) ((b-a)*x + a + b) / 2;
[gp, gw] = grule(n_gp);

phys_jac = rectangle_jacobian_local(ConPts_o);
n_ele_dofs = (pu + 1) * (pv + 1);
nq_elem = n_gp * n_gp;
quad_weights_ref = kron(gw(:), gw(:));

A_value = zeros(NoEs * n_ele_dofs * n_ele_dofs, 1);
row_index = zeros(size(A_value));
column_index = zeros(size(A_value));

global_index = 1;
counter = 1;

% Integrate the sampled potential over each element.
for e = 1:NoEs
    ue = Coordinate(e,1:2);
    ve = Coordinate(e,3:4);
    jac_elem = ((ue(2) - ue(1)) / 2) * ((ve(2) - ve(1)) / 2) * phys_jac;

    row = Element(e,:);

    Nu_vals = zeros(pu + 1, n_gp);
    for i = 1:n_gp
        upar = Fhat(gp(i), ue(1), ue(2));
        Nu_vals(:, i) = bsplinebasis(knotU, pu, upar);
    end

    Nv_vals = zeros(pv + 1, n_gp);
    for j = 1:n_gp
        vpar = Fhat(gp(j), ve(1), ve(2));
        Nv_vals(:, j) = bsplinebasis(knotV, pv, vpar);
    end

    B = zeros(nq_elem, n_ele_dofs);
    qid = 1;
    for i = 1:n_gp
        Nu = Nu_vals(:, i);
        for j = 1:n_gp
            Nv = Nv_vals(:, j).';
            basis_funcs = Nu * Nv;
            B(qid, :) = basis_funcs(:).';
            qid = qid + 1;
        end
    end

    rho_e = rho_samples(counter:counter+nq_elem-1);
    weights_e = rho_e .* (jac_elem * quad_weights_ref);
    Ae = B.' * bsxfun(@times, B, weights_e);
    counter = counter + nq_elem;

    entry_range = (global_index:global_index+n_ele_dofs*n_ele_dofs-1).';
    [rr, cc] = ndgrid(row, row);
    row_index(entry_range) = rr(:);
    column_index(entry_range) = cc(:);
    A_value(entry_range) = Ae(:);
    global_index = global_index + n_ele_dofs * n_ele_dofs;
end

% Assemble and symmetrize the global nonlinear matrix.
N_in = sparse(row_index, column_index, A_value, n_dofs, n_dofs);
N_in = 0.5 * (N_in + N_in');
end

function jac = rectangle_jacobian_local(ConPts)
% Physical Jacobian for an axis-aligned rectangular patch on [0,1]^2.
x = ConPts(:, :, 1);
y = ConPts(:, :, 2);
jac = (max(x(:)) - min(x(:))) * (max(y(:)) - min(y(:)));
end
