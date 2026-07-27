function rho_samples = sample_rho_nurbs(cI, ~, nurbs_refine, n_gp)
% Sample the NURBS density.

Element    = nurbs_refine.Element;
knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
Coordinate = nurbs_refine.Coordinate;
NoEs       = nurbs_refine.NoEs;
pu         = nurbs_refine.pu;
pv         = nurbs_refine.pv;

Fhat = @(x,a,b) ((b-a)*x + a + b) / 2;
[gp, ~] = grule(n_gp);

n_ele_dofs = (pu + 1) * (pv + 1);
nq_elem = n_gp * n_gp;
rho_samples = zeros(NoEs * nq_elem, 1);

counter = 1;
for e = 1:NoEs
    ue = Coordinate(e,1:2);
    ve = Coordinate(e,3:4);

    row = Element(e,:);
    ce = cI(row);

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

    uvals = B * ce;
    rho_samples(counter:counter+nq_elem-1) = real(abs(uvals).^2);
    counter = counter + nq_elem;
end
end
