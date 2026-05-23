function [Vmat, meta] = assemble_NURBS_potential_interp_3D( ...
nurbs_original, nurbs_refine, Vgrid, L, n_gp, opts)
% Assemble the inner IGA potential block from a Cartesian midpoint grid.
arguments
nurbs_original
nurbs_refine
Vgrid
L
n_gp
opts struct
end

use_affine_cube_fast = opts.use_affine_cube_fast;
[supported, support] = build_affine_cube_fast_support_3D(nurbs_original, nurbs_refine, n_gp);

if use_affine_cube_fast && supported
    [Vmat, meta] = assemble_affine_cube_fast_local(support, Vgrid, L);
else
    [Vmat, meta] = assemble_legacy_local(nurbs_original, nurbs_refine, Vgrid, L, n_gp);
end

meta.use_affine_cube_fast_requested = use_affine_cube_fast;
meta.affine_cube_supported = supported;
end

function [Vmat, meta] = assemble_affine_cube_fast_local(support, Vgrid, L)
%Assemble matrices or interface terms for the method.
t_total = tic;

mFFT = size(Vgrid, 1);
dx = L / mFFT;
x1d = -L / 2 + dx / 2 + (0:mFFT-1) * dx;

Element = support.Element;
NoEs = support.NoEs;
n_dofs = support.n_dofs;
uNoEs = support.uNoEs;
vNoEs = support.vNoEs;
wNoEs = support.wNoEs;
n_ele_dofs = support.n_ele_dofs;
n_local_entries = n_ele_dofs * n_ele_dofs;

Vval = zeros(NoEs * n_local_entries, 1);
Iind = zeros(NoEs * n_local_entries, 1);
Jind = zeros(NoEs * n_local_entries, 1);
gidx = 1;

t_interp = 0;
t_local = 0;
t_scatter = 0;

for kw = 1:wNoEs
    wdata = support.w_cache{kw};
    for jv = 1:vNoEs
        vdata = support.v_cache{jv};
        for iu = 1:uNoEs
            e = iu + (jv - 1) * uNoEs + (kw - 1) * uNoEs * vNoEs;
            udata = support.u_cache{iu};
            row = Element(e, :);

            t0 = tic;
            [Phi, wJ, x_pts, y_pts, z_pts] = build_element_tensor_ops_local( ...
                udata, vdata, wdata, support.geom.abs_detDF);
            Vq = interpn(x1d, x1d, x1d, Vgrid, x_pts, y_pts, z_pts, 'linear');
            t_interp = t_interp + toc(t0);

            t0 = tic;
            Phi_v = bsxfun(@times, Phi, reshape(wJ .* Vq, 1, []));
            Vloc = Phi_v * Phi';
            t_local = t_local + toc(t0);

            t0 = tic;
            idx = gidx:(gidx + n_local_entries - 1);
            row_vec = row(:);
            Iind(idx) = repmat(row_vec, n_ele_dofs, 1);
            Jind(idx) = repelem(row_vec, n_ele_dofs, 1);
            Vval(idx) = Vloc(:);
            gidx = gidx + n_local_entries;
            t_scatter = t_scatter + toc(t0);
        end
    end
end

Vmat = sparse(Iind, Jind, Vval, n_dofs, n_dofs);
Vmat = 0.5 * (Vmat + Vmat');

meta = struct();
meta.method = 'affine_cube_fast';
meta.t_interp = t_interp;
meta.t_local = t_local;
meta.t_scatter = t_scatter;
meta.t_total = toc(t_total);
end

function [Phi, wJ, x_pts, y_pts, z_pts] = build_element_tensor_ops_local( ...
udata, vdata, wdata, abs_detDF)

nq_u = numel(udata.phys_pts);
nq_v = numel(vdata.phys_pts);
nq_w = numel(wdata.phys_pts);
nq_total = nq_u * nq_v * nq_w;
n_ele_dofs = size(udata.N, 1) * size(vdata.N, 1) * size(wdata.N, 1);

Phi = zeros(n_ele_dofs, nq_total);
wJ = abs_detDF * kron(wdata.w_line, kron(vdata.w_line, udata.w_line));

[Xu, Xv, Xw] = ndgrid(udata.phys_pts, vdata.phys_pts, wdata.phys_pts);
x_pts = Xu(:);
y_pts = Xv(:);
z_pts = Xw(:);

for kw = 1:nq_w
    for jv = 1:nq_v
        cols = ((kw - 1) * nq_v + (jv - 1)) * nq_u + (1:nq_u);
        Phi(:, cols) = kron(wdata.N(:, kw), kron(vdata.N(:, jv), udata.N));
    end
end
end

function [Vmat, meta] = assemble_legacy_local(nurbs_original, nurbs_refine, Vgrid, L, n_gp)
%Assemble matrices or interface terms for the method.
t_total = tic;

ConPts_o = nurbs_original.ConPts;
weights_o = nurbs_original.weights;
knotU_o = nurbs_original.knotU;
knotV_o = nurbs_original.knotV;
knotW_o = nurbs_original.knotW;
pu_o = nurbs_original.pu;
pv_o = nurbs_original.pv;
pw_o = nurbs_original.pw;

Element = nurbs_refine.Element;
Coordinate = nurbs_refine.Coordinate;
NoEs = nurbs_refine.NoEs;
n_dofs = nurbs_refine.n_dofs;

knotU = nurbs_refine.Ubar;
knotV = nurbs_refine.Vbar;
knotW = nurbs_refine.Wbar;
pu = nurbs_refine.pu;
pv = nurbs_refine.pv;
pw = nurbs_refine.pw;

mFFT = size(Vgrid, 1);
dx = L / mFFT;
x1d = -L / 2 + dx / 2 + (0:mFFT-1) * dx;

[gp, gw] = grule(n_gp);
Fhat = @(x, a, b) ((b - a) * x + a + b) / 2;

n_ele_dofs = (pu + 1) * (pv + 1) * (pw + 1);
Vval = zeros(NoEs * n_ele_dofs * n_ele_dofs, 1);
Iind = Vval;
Jind = Vval;
gidx = 1;

for e = 1:NoEs
    Ve = zeros(n_ele_dofs, n_ele_dofs);

    ue = Coordinate(e, 1:2);
    ve = Coordinate(e, 3:4);
    we = Coordinate(e, 5:6);
    J1 = (ue(2) - ue(1)) / 2;
    J2 = (ve(2) - ve(1)) / 2;
    J3 = (we(2) - we(1)) / 2;
    row = Element(e, :);

    for i = 1:n_gp
        u = Fhat(gp(i), ue(1), ue(2));
        Uders = bspbasisDers(knotU, pu, u, 1);
        Nu = Uders(1, :)';

        for j = 1:n_gp
            v = Fhat(gp(j), ve(1), ve(2));
            Vders = bspbasisDers(knotV, pv, v, 1);
            Nv = Vders(1, :)';

            for k = 1:n_gp
                w = Fhat(gp(k), we(1), we(2));
                Wders = bspbasisDers(knotW, pw, w, 1);
                Nw = Wders(1, :)';

                [F, DF] = NurbsVolume(ConPts_o, weights_o, ...
                    knotU_o, pu_o, u, knotV_o, pv_o, v, knotW_o, pw_o, w);

                basis = zeros(n_ele_dofs, 1);
                cnt = 0;
                for kk = 1:(pw + 1)
                    for jj = 1:(pv + 1)
                        for ii = 1:(pu + 1)
                            cnt = cnt + 1;
                            basis(cnt) = Nu(ii) * Nv(jj) * Nw(kk);
                        end
                    end
                end

                Vq = interpn(x1d, x1d, x1d, Vgrid, F(1), F(2), F(3), 'linear');
                Jac = J1 * J2 * J3 * abs(det(DF)) * gw(i) * gw(j) * gw(k);

                Ve = Ve + Vq * (basis * basis') * Jac;
            end
        end
    end

    for i1 = 1:n_ele_dofs
        for j1 = 1:n_ele_dofs
            Iind(gidx) = row(i1);
            Jind(gidx) = row(j1);
            Vval(gidx) = Ve(i1, j1);
            gidx = gidx + 1;
        end
    end
end

Vmat = sparse(Iind, Jind, Vval, n_dofs, n_dofs);
Vmat = 0.5 * (Vmat + Vmat');

meta = struct();
meta.method = 'legacy_general';
meta.t_total = toc(t_total);
end
