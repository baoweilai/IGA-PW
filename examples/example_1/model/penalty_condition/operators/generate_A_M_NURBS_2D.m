function [A, M, meta] = generate_A_M_NURBS_2D( ...
nurbs_original, nurbs_refine, p_Vr, n_pw_Vr, L, n_gp, Example, varargin)

opts = struct();
if ~isempty(varargin)
    opts = varargin{1};
end

use_square_nurbs_fast = opts.use_square_nurbs_fast;
[square_supported, geom] = detect_affine_square_geometry_local(nurbs_original);

if use_square_nurbs_fast && square_supported
    [A, M, meta] = generate_A_M_NURBS_2D_square_fast_local( ...
        nurbs_refine, p_Vr, n_pw_Vr, L, n_gp, geom);
else
    [A, M, meta] = generate_A_M_NURBS_2D_legacy_local( ...
        nurbs_original, nurbs_refine, p_Vr, n_pw_Vr, L, n_gp, Example);
end

meta.use_square_nurbs_fast_requested = use_square_nurbs_fast;
meta.square_fast_supported = square_supported;
if square_supported
    meta.square_geometry = geom;
end
end

function [A, M, meta] = generate_A_M_NURBS_2D_square_fast_local( ...
nurbs_refine, p_Vr, n_pw_Vr, L, n_gp, geom)

t_total = tic;
t_precompute = tic;

Element = nurbs_refine.Element;
knotU = nurbs_refine.Ubar;
knotV = nurbs_refine.Vbar;
Coordinate = nurbs_refine.Coordinate;
UBreaks = nurbs_refine.UBreaks;
VBreaks = nurbs_refine.VBreaks;
NoEs = nurbs_refine.NoEs;
n_dofs = nurbs_refine.n_dofs;
uNoEs = nurbs_refine.uNoEs;
vNoEs = nurbs_refine.vNoEs;

pu = nurbs_refine.pu;
pv = nurbs_refine.pv;
[gp, gw] = grule(n_gp);
n_ele_dofs = (pu + 1) * (pv + 1);
n_local_entries = n_ele_dofs * n_ele_dofs;

u_cache = precompute_basis_line_local(UBreaks, knotU, pu, gp, gw, geom.x0, geom.sx);
v_cache = precompute_basis_line_local(VBreaks, knotV, pv, gp, gw, geom.y0, geom.sy);

t_precompute = toc(t_precompute);

A_value = zeros(NoEs * n_local_entries, 1);
M_value = zeros(NoEs * n_local_entries, 1);
row_index = zeros(NoEs * n_local_entries, 1);
column_index = zeros(NoEs * n_local_entries, 1);

abs_detDF = abs(geom.sx * geom.sy);
global_index = 1;
t_potential = 0;
t_local_matrix = 0;
t_scatter = 0;

for jv = 1:vNoEs
    vdata = v_cache{jv};
    for iu = 1:uNoEs
        e = iu + (jv - 1) * uNoEs;
        udata = u_cache{iu};
        row = Element(e, :);

        t0 = tic;
        [Phi, Gx, Gy, wJ, x_pts, y_pts] = build_element_tensor_operators_local( ...
            udata, vdata, abs_detDF);
        Vvals = Vr_2D_Example_1_batch_local(p_Vr, n_pw_Vr, L, x_pts, y_pts);
        t_potential = t_potential + toc(t0);

        t0 = tic;
        wJ_row = reshape(wJ, 1, []);
        Phi_w = bsxfun(@times, Phi, wJ_row);
        Mloc = Phi_w * Phi';

        Gx_w = bsxfun(@times, Gx, wJ_row);
        Gy_w = bsxfun(@times, Gy, wJ_row);
        Kloc = 0.5 * (Gx_w * Gx' + Gy_w * Gy');

        Phi_v = bsxfun(@times, Phi, reshape(wJ .* Vvals, 1, []));
        Aloc = Kloc + Phi_v * Phi';
        t_local_matrix = t_local_matrix + toc(t0);

        t0 = tic;
        idx = global_index:(global_index + n_local_entries - 1);
        row_vec = row(:);
        row_index(idx) = repmat(row_vec, n_ele_dofs, 1);
        column_index(idx) = repelem(row_vec, n_ele_dofs, 1);
        A_value(idx) = Aloc(:);
        M_value(idx) = Mloc(:);
        global_index = global_index + n_local_entries;
        t_scatter = t_scatter + toc(t0);
    end
end

A = sparse(row_index, column_index, A_value, n_dofs, n_dofs);
M = sparse(row_index, column_index, M_value, n_dofs, n_dofs);

meta = struct();
meta.method = 'square_fast';
meta.n_gp = n_gp;
meta.t_precompute = t_precompute;
meta.t_potential = t_potential;
meta.t_local_matrix = t_local_matrix;
meta.t_scatter = t_scatter;
meta.t_total = toc(t_total);
meta.geometry = geom;
meta.n_elements = size(Coordinate, 1);
end

function [Phi, Gx, Gy, wJ, x_pts, y_pts] = build_element_tensor_operators_local( ...
udata, vdata, abs_detDF)

n_gp_u = numel(udata.phys_pts);
n_gp_v = numel(vdata.phys_pts);
n_gp_total = n_gp_u * n_gp_v;
n_ele_dofs = size(udata.N, 1) * size(vdata.N, 1);

Phi = zeros(n_ele_dofs, n_gp_total);
Gx = zeros(n_ele_dofs, n_gp_total);
Gy = zeros(n_ele_dofs, n_gp_total);

x_pts = repmat(udata.phys_pts(:), n_gp_v, 1);
y_pts = kron(vdata.phys_pts(:), ones(n_gp_u, 1));
wJ = abs_detDF * kron(vdata.w_line, udata.w_line);

for j = 1:n_gp_v
    cols = (j - 1) * n_gp_u + (1:n_gp_u);
    Phi(:, cols) = kron(vdata.N(:, j), udata.N);
    Gx(:, cols) = kron(vdata.N(:, j), udata.D_phys);
    Gy(:, cols) = kron(vdata.D_phys(:, j), udata.N);
end
end

function cache = precompute_basis_line_local(Breaks, knot, p, gp, gw, phys0, physScale)
%Compute basis line.

n_ele = numel(Breaks) - 1;
n_gp = numel(gp);
cache = cell(n_ele, 1);

for e = 1:n_ele
    interval = Breaks(e:e + 1);
    J = (interval(2) - interval(1)) / 2;
    param_pts = ((interval(2) - interval(1)) * gp + interval(1) + interval(2)) / 2;

    N = zeros(p + 1, n_gp);
    D_phys = zeros(p + 1, n_gp);
    for ig = 1:n_gp
        ders = bspbasisDers(knot, p, param_pts(ig), 1);
        N(:, ig) = ders(1, :)';
        D_phys(:, ig) = ders(2, :)' / physScale;
    end

    entry = struct();
    entry.N = N;
    entry.D_phys = D_phys;
    entry.w_line = J * gw(:);
    entry.phys_pts = phys0 + physScale * param_pts(:);
    cache{e} = entry;
end
end

function Vr = Vr_2D_Example_1_batch_local(p, n_pw, L, x, y)
%Compute 2D example 1 batch.

p = p * (2 * pi / L);
pts = [x(:), y(:)];
r_norm = hypot(pts(:, 1), pts(:, 2));
alpha = 5;

p_norm = hypot(p(:, 1), p(:, 2));
mask = p_norm > 0;
coeff = erfc(p_norm(mask) / (2 * alpha)) ./ p_norm(mask);

if any(mask)
    phase = pts * p(mask, :)';
    s = exp(1i * phase) * coeff;
else
    s = zeros(size(pts, 1), 1);
end

Vr = -erfc(alpha * r_norm) ./ r_norm - 2 * pi * s / (L * L) + 2 * alpha / sqrt(pi);
Vr = reshape(Vr, [], 1);

if numel(Vr) ~= numel(x) || n_pw ~= size(p, 1)
    error('Vr_2D_Example_1_batch_local: inconsistent input dimensions.');
end
end

function [A, M, meta] = generate_A_M_NURBS_2D_legacy_local( ...
nurbs_original, nurbs_refine, p_Vr, n_pw_Vr, L, n_gp, Example) %#ok<INUSD>

t_total = tic;

ConPts_o = nurbs_original.ConPts;
weights_o = nurbs_original.weights;
knotU_o = nurbs_original.knotU;
knotV_o = nurbs_original.knotV;
pu_o = nurbs_original.pu;
pv_o = nurbs_original.pv;
Element = nurbs_refine.Element;
knotU = nurbs_refine.Ubar;
knotV = nurbs_refine.Vbar;
Coordinate = nurbs_refine.Coordinate;
Fhat = @(x, a, b) ((b - a) * x + a + b) / 2;
NoEs = nurbs_refine.NoEs;
n_dofs = nurbs_refine.n_dofs;

pu = nurbs_refine.pu;
pv = nurbs_refine.pv;
[gp, gw] = grule(n_gp);
n_ele_dofs = (pu + 1) * (pv + 1);

Ae = zeros(n_ele_dofs, n_ele_dofs);
Me = zeros(n_ele_dofs, n_ele_dofs);

A_value = zeros(NoEs * n_ele_dofs * n_ele_dofs, 1);
M_value = zeros(NoEs * n_ele_dofs * n_ele_dofs, 1);
row_index = A_value;
column_index = A_value;
global_index = 1;

for e = 1:NoEs
    Ae = 0 * Ae;
    Me = 0 * Me;
    ue = Coordinate(e, 1:2);
    J1 = (ue(2) - ue(1)) / 2;
    ve = Coordinate(e, 3:4);
    J2 = (ve(2) - ve(1)) / 2;
    row = Element(e, :);
    for i = 1:n_gp
        u = Fhat(gp(i), ue(1), ue(2));
        uders = bspbasisDers(knotU, pu, u, 1);
        Nu = uders(1, :)';
        DNu = uders(2, :)';
        for j = 1:n_gp
            v = Fhat(gp(j), ve(1), ve(2));
            vders = bspbasisDers(knotV, pv, v, 1);
            Nv = vders(1, :);
            DNv = vders(2, :);
            [F, DF] = NurbsSurface(ConPts_o, weights_o, knotU_o, pu_o, u, knotV_o, pv_o, v);
            basis_funcs = Nu * Nv;
            basis_funcs = basis_funcs(:);
            DNu_v = DNu * Nv;
            DNu_v = DNu_v(:);
            DNv_u = Nu * DNv;
            DNv_u = DNv_u(:);
            basis_grad = [DNu_v, DNv_u] / DF;
            Jacobian = J1 * J2 * abs(det(DF)) * gw(i) * gw(j);
            Mass_mat_ele = basis_funcs * basis_funcs' * Jacobian;
            V_ext = Vr_2D_Example_1(p_Vr, L, n_pw_Vr, F(1), F(2));
            Ae = Ae + basis_grad * basis_grad' * Jacobian / 2 + V_ext * Mass_mat_ele;
            Me = Me + Mass_mat_ele;
        end
    end

    for i1 = 1:n_ele_dofs
        for j1 = 1:n_ele_dofs
            row_index(global_index) = row(i1);
            column_index(global_index) = row(j1);
            A_value(global_index) = Ae(i1, j1);
            M_value(global_index) = Me(i1, j1);
            global_index = global_index + 1;
        end
    end
end

A = sparse(row_index, column_index, A_value, n_dofs, n_dofs);
M = sparse(row_index, column_index, M_value, n_dofs, n_dofs);

meta = struct();
meta.method = 'legacy_general';
meta.t_total = toc(t_total);
end

function [tf, geom] = detect_affine_square_geometry_local(nurbs_original)
%Detect affine geometry for a square patch.

geom = struct('x0', 0, 'y0', 0, 'sx', 0, 'sy', 0);
tf = false;

ConPts = nurbs_original.ConPts;
weights = nurbs_original.weights;

[n_u, n_v, dim] = size(ConPts);
if dim ~= 2 || n_u < 2 || n_v < 2 || any(abs(weights(:) - 1) > 1e-14)
    return;
end

x_grid = ConPts(:, :, 1);
y_grid = ConPts(:, :, 2);

x_vary_v = max(abs(diff(x_grid, 1, 2)), [], 1);
x_vary_v = max(x_vary_v(:));
y_vary_u = max(abs(diff(y_grid, 1, 1)), [], 2);
y_vary_u = max(y_vary_u(:));
x_vary_u = max(abs(diff(x_grid, 1, 1)), [], 2);
x_vary_u = max(x_vary_u(:));
y_vary_v = max(abs(diff(y_grid, 1, 2)), [], 1);
y_vary_v = max(y_vary_v(:));

if x_vary_v > 1e-12 || y_vary_u > 1e-12
    return;
end
if x_vary_u <= 1e-14 || y_vary_v <= 1e-14
    return;
end

geom.x0 = x_grid(1, 1);
geom.y0 = y_grid(1, 1);
geom.sx = x_grid(end, 1) - x_grid(1, 1);
geom.sy = y_grid(1, end) - y_grid(1, 1);
geom.detDF = geom.sx * geom.sy;

tf = true;
end
