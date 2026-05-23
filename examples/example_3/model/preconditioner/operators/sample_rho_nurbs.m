function rho_samples = sample_rho_nurbs(cI, ~, nurbs_refine, n_gp)
%Sample the NURBS density.

Element    = nurbs_refine.Element;
knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
Coordinate = nurbs_refine.Coordinate;
NoEs       = nurbs_refine.NoEs;
pu         = nurbs_refine.pu;
pv         = nurbs_refine.pv;

Fhat = @(x,a,b) ((b-a)*x + a + b) / 2;
[gp, ~] = grule(n_gp);

rho_samples = zeros(NoEs * n_gp * n_gp, 1);

counter = 1;
for e = 1:NoEs
    ue = Coordinate(e,1:2);
    ve = Coordinate(e,3:4);

    row = Element(e,:);
    ce = cI(row);

    for i = 1:n_gp
        upar  = Fhat(gp(i), ue(1), ue(2));
        uders = bspbasisDers(knotU, pu, upar, 1);
        Nu    = uders(1,:).';

        for j = 1:n_gp
            vpar  = Fhat(gp(j), ve(1), ve(2));
            vders = bspbasisDers(knotV, pv, vpar, 1);
            Nv    = vders(1,:);

            basis_funcs = Nu * Nv;
            basis_funcs = basis_funcs(:);

            uval = basis_funcs.' * ce;
            rho_samples(counter) = real(abs(uval)^2);
            counter = counter + 1;
        end
    end
end
end
