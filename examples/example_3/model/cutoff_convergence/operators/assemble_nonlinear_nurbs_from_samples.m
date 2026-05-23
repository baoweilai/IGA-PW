function N_in = assemble_nonlinear_nurbs_from_samples( ...
nurbs_original, nurbs_refine, rho_samples, n_gp)

ConPts_o   = nurbs_original.ConPts;
weights_o  = nurbs_original.weights;
knotU_o    = nurbs_original.knotU;
knotV_o    = nurbs_original.knotV;
pu_o       = nurbs_original.pu;
pv_o       = nurbs_original.pv;

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
    error('assemble_nonlinear_nurbs_from_samples: rho_samples length mismatch.');
end

Fhat = @(x,a,b) ((b-a)*x + a + b) / 2;
[gp, gw] = grule(n_gp);

n_ele_dofs = (pu + 1) * (pv + 1);
Ae = zeros(n_ele_dofs, n_ele_dofs);

A_value = zeros(NoEs * n_ele_dofs * n_ele_dofs, 1);
row_index = zeros(size(A_value));
column_index = zeros(size(A_value));

global_index = 1;
counter = 1;

for e = 1:NoEs
    Ae(:) = 0;

    ue = Coordinate(e,1:2);
    ve = Coordinate(e,3:4);
    J1 = (ue(2) - ue(1)) / 2;
    J2 = (ve(2) - ve(1)) / 2;

    row = Element(e,:);

    for i = 1:n_gp
        upar  = Fhat(gp(i), ue(1), ue(2));
        uders = bspbasisDers(knotU, pu, upar, 1);
        Nu    = uders(1,:).';

        for j = 1:n_gp
            vpar  = Fhat(gp(j), ve(1), ve(2));
            vders = bspbasisDers(knotV, pv, vpar, 1);
            Nv    = vders(1,:);

            [~, DF] = NurbsSurface(ConPts_o, weights_o, knotU_o, pu_o, upar, ...
                knotV_o, pv_o, vpar);

            basis_funcs = Nu * Nv;
            basis_funcs = basis_funcs(:);

            Jacobian = J1 * J2 * abs(det(DF)) * gw(i) * gw(j);

            Ae = Ae + rho_samples(counter) * (basis_funcs * basis_funcs.') * Jacobian;
            counter = counter + 1;
        end
    end

    for i1 = 1:n_ele_dofs
        for j1 = 1:n_ele_dofs
            row_index(global_index) = row(i1);
            column_index(global_index) = row(j1);
            A_value(global_index) = Ae(i1,j1);
            global_index = global_index + 1;
        end
    end
end

N_in = sparse(row_index, column_index, A_value, n_dofs, n_dofs);
N_in = 0.5 * (N_in + N_in');
end
