function [K, M, meta] = generate_K_M_NURBS_3D(nurbs_original, nurbs_refine, n_gp, opts)
% Assemble 3-D NURBS stiffness and mass matrices.
arguments
    nurbs_original
    nurbs_refine
    n_gp
    opts struct
end

% Select affine tensor assembly when the geometry supports it.
use_affine_cube_fast = opts.use_affine_cube_fast;
[supported, support] = build_affine_cube_fast_support_3D(nurbs_original, nurbs_refine, n_gp);

if use_affine_cube_fast && supported
    [K, M, meta] = generate_K_M_affine_cube_fast_local(support);
else
    [K, M, meta] = generate_K_M_legacy_local(nurbs_original, nurbs_refine, n_gp);
end

meta.use_affine_cube_fast_requested = use_affine_cube_fast;
meta.affine_cube_supported = supported;
end

function [K, M, meta] = generate_K_M_affine_cube_fast_local(support)
% Assemble stiffness and mass matrices with affine tensor products.
t_total = tic;

% Allocate global triplets and timing accumulators.
Element = support.Element;
NoEs = support.NoEs;
n_dofs = support.n_dofs;
uNoEs = support.uNoEs;
vNoEs = support.vNoEs;
wNoEs = support.wNoEs;
n_ele_dofs = support.n_ele_dofs;
n_local_entries = n_ele_dofs * n_ele_dofs;

Kval = zeros(NoEs * n_local_entries, 1);
Mval = zeros(NoEs * n_local_entries, 1);
Iind = zeros(NoEs * n_local_entries, 1);
Jind = zeros(NoEs * n_local_entries, 1);

gidx = 1;
t_local = 0;
t_scatter = 0;

% Build and scatter each tensor-product element block.
for kw = 1:wNoEs
    wdata = support.w_cache{kw};
    for jv = 1:vNoEs
        vdata = support.v_cache{jv};
        for iu = 1:uNoEs
            e = iu + (jv - 1) * uNoEs + (kw - 1) * uNoEs * vNoEs;
            udata = support.u_cache{iu};
            row = Element(e, :);

            t0 = tic;
            [Phi, Gx, Gy, Gz, wJ] = build_element_tensor_operators_local( ...
                udata, vdata, wdata, support.geom.abs_detDF);

            wJ_row = reshape(wJ, 1, []);
            Phi_w = bsxfun(@times, Phi, wJ_row);
            Mloc = Phi_w * Phi';

            Gx_w = bsxfun(@times, Gx, wJ_row);
            Gy_w = bsxfun(@times, Gy, wJ_row);
            Gz_w = bsxfun(@times, Gz, wJ_row);
            Kloc = 0.5 * (Gx_w * Gx' + Gy_w * Gy' + Gz_w * Gz');
            t_local = t_local + toc(t0);

            t0 = tic;
            idx = gidx:(gidx + n_local_entries - 1);
            row_vec = row(:);
            Iind(idx) = repmat(row_vec, n_ele_dofs, 1);
            Jind(idx) = repelem(row_vec, n_ele_dofs, 1);
            Kval(idx) = Kloc(:);
            Mval(idx) = Mloc(:);
            gidx = gidx + n_local_entries;
            t_scatter = t_scatter + toc(t0);
        end
    end
end

% Assemble the global sparse matrices and timing data.
K = sparse(Iind, Jind, Kval, n_dofs, n_dofs);
M = sparse(Iind, Jind, Mval, n_dofs, n_dofs);

meta = struct();
meta.method = 'affine_cube_fast';
meta.t_local = t_local;
meta.t_scatter = t_scatter;
meta.t_total = toc(t_total);
end

function [Phi, Gx, Gy, Gz, wJ] = build_element_tensor_operators_local( ...
udata, vdata, wdata, abs_detDF)

% Build tensor-product basis, gradient, and quadrature operators for one element.
nq_u = numel(udata.phys_pts);
nq_v = numel(vdata.phys_pts);
nq_w = numel(wdata.phys_pts);
nq_total = nq_u * nq_v * nq_w;
n_ele_dofs = size(udata.N, 1) * size(vdata.N, 1) * size(wdata.N, 1);

Phi = zeros(n_ele_dofs, nq_total);
Gx = zeros(n_ele_dofs, nq_total);
Gy = zeros(n_ele_dofs, nq_total);
Gz = zeros(n_ele_dofs, nq_total);

wJ = abs_detDF * kron(wdata.w_line, kron(vdata.w_line, udata.w_line));

for kw = 1:nq_w
    for jv = 1:nq_v
        cols = ((kw - 1) * nq_v + (jv - 1)) * nq_u + (1:nq_u);
        Phi(:, cols) = kron(wdata.N(:, kw), kron(vdata.N(:, jv), udata.N));
        Gx(:, cols) = kron(wdata.N(:, kw), kron(vdata.N(:, jv), udata.D_phys));
        Gy(:, cols) = kron(wdata.N(:, kw), kron(vdata.D_phys(:, jv), udata.N));
        Gz(:, cols) = kron(wdata.D_phys(:, kw), kron(vdata.N(:, jv), udata.N));
    end
end
end

function [K, M, meta] = generate_K_M_legacy_local(nurbs_original, nurbs_refine, n_gp)
% Assemble stiffness and mass matrices by element quadrature.
t_total = tic;

% Read the geometry and refined-mesh data.
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

% Prepare quadrature and sparse matrix triplets.
[gp, gw] = grule(n_gp);
Fhat = @(x, a, b) ((b - a) * x + a + b) / 2;

n_ele_dofs = (pu + 1) * (pv + 1) * (pw + 1);
Kval = zeros(NoEs * n_ele_dofs * n_ele_dofs, 1);
Mval = zeros(NoEs * n_ele_dofs * n_ele_dofs, 1);
Iind = Kval;
Jind = Kval;

gidx = 1;

% Integrate each element in physical coordinates.
for e = 1:NoEs
    Ke = zeros(n_ele_dofs, n_ele_dofs);
    Me = zeros(n_ele_dofs, n_ele_dofs);

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
        DNu = Uders(2, :)';

        for j = 1:n_gp
            v = Fhat(gp(j), ve(1), ve(2));
            Vders = bspbasisDers(knotV, pv, v, 1);
            Nv = Vders(1, :)';
            DNv = Vders(2, :)';

            for k = 1:n_gp
                w = Fhat(gp(k), we(1), we(2));
                Wders = bspbasisDers(knotW, pw, w, 1);
                Nw = Wders(1, :)';
                DNw = Wders(2, :)';

                [~, DF] = NurbsVolume(ConPts_o, weights_o, ...
                    knotU_o, pu_o, u, knotV_o, pv_o, v, knotW_o, pw_o, w);

                basis = zeros(n_ele_dofs, 1);
                dU = basis;
                dV = basis;
                dW = basis;

                cnt = 0;
                for kk = 1:(pw + 1)
                    for jj = 1:(pv + 1)
                        for ii = 1:(pu + 1)
                            cnt = cnt + 1;
                            basis(cnt) = Nu(ii) * Nv(jj) * Nw(kk);
                            dU(cnt) = DNu(ii) * Nv(jj) * Nw(kk);
                            dV(cnt) = Nu(ii) * DNv(jj) * Nw(kk);
                            dW(cnt) = Nu(ii) * Nv(jj) * DNw(kk);
                        end
                    end
                end

                grad = [dU, dV, dW] / DF;
                Jac = J1 * J2 * J3 * abs(det(DF)) * gw(i) * gw(j) * gw(k);

                Me = Me + (basis * basis') * Jac;
                Ke = Ke + (grad * grad') * (Jac / 2);
            end
        end
    end

    % Store the element blocks in global triplets.
    for i1 = 1:n_ele_dofs
        for j1 = 1:n_ele_dofs
            Iind(gidx) = row(i1);
            Jind(gidx) = row(j1);
            Kval(gidx) = Ke(i1, j1);
            Mval(gidx) = Me(i1, j1);
            gidx = gidx + 1;
        end
    end
end

% Assemble the global sparse matrices and timing data.
K = sparse(Iind, Jind, Kval, n_dofs, n_dofs);
M = sparse(Iind, Jind, Mval, n_dofs, n_dofs);

meta = struct();
meta.method = 'legacy_general';
meta.t_total = toc(t_total);
end
