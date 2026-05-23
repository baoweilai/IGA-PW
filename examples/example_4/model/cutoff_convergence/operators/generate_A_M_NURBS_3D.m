function [A,M] = generate_A_M_NURBS_3D(nurbs_original, nurbs_refine, p_Vr, n_pw_Vr, L, n_gp, Example)
%Generate a m NURBS 3D.

ConPts_o   = nurbs_original.ConPts;
weights_o  = nurbs_original.weights;
knotU_o    = nurbs_original.knotU;
knotV_o    = nurbs_original.knotV;
knotW_o    = nurbs_original.knotW;
pu_o       = nurbs_original.pu;
pv_o       = nurbs_original.pv;
pw_o       = nurbs_original.pw;

Element    = nurbs_refine.Element;
knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
knotW      = nurbs_refine.Wbar;
Coordinate = nurbs_refine.Coordinate;
NoEs       = nurbs_refine.NoEs;
n_dofs     = nurbs_refine.n_dofs;

pu = nurbs_refine.pu;
pv = nurbs_refine.pv;
pw = nurbs_refine.pw;

[gp,gw] = grule(n_gp);
Fhat = @(x,a,b) ((b-a)*x + a + b)/2;

n_ele_dofs = (pu+1)*(pv+1)*(pw+1);

A_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
M_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
row_idx = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
col_idx = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);

gidx = 1;

for e = 1:NoEs
    Ae = zeros(n_ele_dofs,n_ele_dofs);
    Me = zeros(n_ele_dofs,n_ele_dofs);

    ue = Coordinate(e,1:2);
    ve = Coordinate(e,3:4);
    we = Coordinate(e,5:6);

    J1 = (ue(2) - ue(1))/2;
    J2 = (ve(2) - ve(1))/2;
    J3 = (we(2) - we(1))/2;

    row = Element(e,:);

    for i = 1:n_gp
        u = Fhat(gp(i), ue(1), ue(2));
        Uders = bspbasisDers(knotU, pu, u, 1);
        Nu  = Uders(1,:)';
        DNu = Uders(2,:)';

        for j = 1:n_gp
            v = Fhat(gp(j), ve(1), ve(2));
            Vders = bspbasisDers(knotV, pv, v, 1);
            Nv  = Vders(1,:)';
            DNv = Vders(2,:)';

            for k = 1:n_gp
                w = Fhat(gp(k), we(1), we(2));
                Wders = bspbasisDers(knotW, pw, w, 1);
                Nw  = Wders(1,:)';
                DNw = Wders(2,:)';

                [F, DF] = NurbsVolume(ConPts_o, weights_o, ...
                    knotU_o, pu_o, u, knotV_o, pv_o, v, knotW_o, pw_o, w);

                basis = zeros(n_ele_dofs,1);
                dU = basis;
                dV = basis;
                dW = basis;

                cnt = 0;
                for kk = 1:(pw+1)
                    for jj = 1:(pv+1)
                        for ii = 1:(pu+1)
                            cnt = cnt + 1;
                            basis(cnt) = Nu(ii)  * Nv(jj)  * Nw(kk);
                            dU(cnt)    = DNu(ii) * Nv(jj)  * Nw(kk);
                            dV(cnt)    = Nu(ii)  * DNv(jj) * Nw(kk);
                            dW(cnt)    = Nu(ii)  * Nv(jj)  * DNw(kk);
                        end
                    end
                end

                grad = [dU, dV, dW] / DF;
                Jac = J1 * J2 * J3 * abs(det(DF)) * gw(i) * gw(j) * gw(k);

                Mass_ele = (basis * basis') * Jac;

                if strcmp(Example, 'Example_1_3D') || strcmp(Example, 'Example_1')
                    V_ext = Vr_3D_Example_1(p_Vr, L, n_pw_Vr, F(1), F(2), F(3));
                else
                    error('generate_A_M_NURBS_3D currently only supports Example_1 / Example_1_3D.');
                end

                Ae = Ae + (grad * grad') * (Jac/2) + V_ext * Mass_ele;
                Me = Me + Mass_ele;
            end
        end
    end

    for i1 = 1:n_ele_dofs
        for j1 = 1:n_ele_dofs
            row_idx(gidx) = row(i1);
            col_idx(gidx) = row(j1);
            A_value(gidx) = Ae(i1,j1);
            M_value(gidx) = Me(i1,j1);
            gidx = gidx + 1;
        end
    end
end

A = sparse(row_idx, col_idx, A_value, n_dofs, n_dofs);
M = sparse(row_idx, col_idx, M_value, n_dofs, n_dofs);

A = 0.5 * (A + A');
M = 0.5 * (M + M');
end
