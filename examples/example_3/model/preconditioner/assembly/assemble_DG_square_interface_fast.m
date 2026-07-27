function [P, S, meta] = assemble_DG_square_interface_fast( ...
nurbs_refine, pw_index, plane_wave_dofs_index, L, a, n_dofs)
% Fast DG assembly for the axis-aligned square inner patch [-a,a]^2.
% This exploits the straight-edge geometry:
%   1) the physical map is affine, so normals/Jacobians are constant;
%   2) the PW-PW edge blocks can be assembled once per whole edge.

t_total = tic;

alpha = 2 * pi / L;
Omega_area = L * L;
pwx = pw_index(:, 1);
pwy = pw_index(:, 2);
n_pw = size(pw_index, 1);
pw_dofs = plane_wave_dofs_index(:);

P = sparse(n_dofs, n_dofs);
S = sparse(n_dofs, n_dofs);
Ppw = complex(zeros(n_pw, n_pw));
Spw = complex(zeros(n_pw, n_pw));

diffx = pwx - pwx.';
diffy = pwy - pwy.';
Kx = interval_ft_general_local(diffx, -a, a, alpha);
Ky = interval_ft_general_local(diffy, -a, a, alpha);

inv_df = 1 / (2 * a);

t_edges = tic;

% Bottom edge: v = 0, y = -a, n = [0; -1]
Vders_bottom = bspbasisDers(nurbs_refine.Vbar, nurbs_refine.pv, 0, 1);
[Pedge, Sedge, Pmm, Smm, edgeMetaBottom] = assemble_horizontal_edge_local( ...
    nurbs_refine.UBreaks, nurbs_refine.Ubar, nurbs_refine.pu, ...
    Vders_bottom(1, 1:2), Vders_bottom(2, 1:2), ...
    nurbs_refine.bottom_dofs_2_layers, -a, -1, ...
    pwx, pwy, alpha, Omega_area, a, inv_df, Kx);
P(nurbs_refine.bottom_dofs_2_layers, nurbs_refine.bottom_dofs_2_layers) = ...
    P(nurbs_refine.bottom_dofs_2_layers, nurbs_refine.bottom_dofs_2_layers) + Pedge.Ppp;
P(nurbs_refine.bottom_dofs_2_layers, pw_dofs) = ...
    P(nurbs_refine.bottom_dofs_2_layers, pw_dofs) + Pedge.Ppm;
P(pw_dofs, nurbs_refine.bottom_dofs_2_layers) = ...
    P(pw_dofs, nurbs_refine.bottom_dofs_2_layers) + Pedge.Pmp;
S(nurbs_refine.bottom_dofs_2_layers, nurbs_refine.bottom_dofs_2_layers) = ...
    S(nurbs_refine.bottom_dofs_2_layers, nurbs_refine.bottom_dofs_2_layers) + Sedge.Spp;
S(nurbs_refine.bottom_dofs_2_layers, pw_dofs) = ...
    S(nurbs_refine.bottom_dofs_2_layers, pw_dofs) + Sedge.Spm;
S(pw_dofs, nurbs_refine.bottom_dofs_2_layers) = ...
    S(pw_dofs, nurbs_refine.bottom_dofs_2_layers) + Sedge.Smp;
Ppw = Ppw + Pmm;
Spw = Spw + Smm;

% Top edge: v = 1, y = a, n = [0; 1]
Vders_top = bspbasisDers(nurbs_refine.Vbar, nurbs_refine.pv, 1, 1);
[Pedge, Sedge, Pmm, Smm, edgeMetaTop] = assemble_horizontal_edge_local( ...
    nurbs_refine.UBreaks, nurbs_refine.Ubar, nurbs_refine.pu, ...
    Vders_top(1, end-1:end), Vders_top(2, end-1:end), ...
    nurbs_refine.top_dofs_2_layers, a, 1, ...
    pwx, pwy, alpha, Omega_area, a, inv_df, Kx);
P(nurbs_refine.top_dofs_2_layers, nurbs_refine.top_dofs_2_layers) = ...
    P(nurbs_refine.top_dofs_2_layers, nurbs_refine.top_dofs_2_layers) + Pedge.Ppp;
P(nurbs_refine.top_dofs_2_layers, pw_dofs) = ...
    P(nurbs_refine.top_dofs_2_layers, pw_dofs) + Pedge.Ppm;
P(pw_dofs, nurbs_refine.top_dofs_2_layers) = ...
    P(pw_dofs, nurbs_refine.top_dofs_2_layers) + Pedge.Pmp;
S(nurbs_refine.top_dofs_2_layers, nurbs_refine.top_dofs_2_layers) = ...
    S(nurbs_refine.top_dofs_2_layers, nurbs_refine.top_dofs_2_layers) + Sedge.Spp;
S(nurbs_refine.top_dofs_2_layers, pw_dofs) = ...
    S(nurbs_refine.top_dofs_2_layers, pw_dofs) + Sedge.Spm;
S(pw_dofs, nurbs_refine.top_dofs_2_layers) = ...
    S(pw_dofs, nurbs_refine.top_dofs_2_layers) + Sedge.Smp;
Ppw = Ppw + Pmm;
Spw = Spw + Smm;

% Left edge: u = 0, x = -a, n = [-1; 0]
Uders_left = bspbasisDers(nurbs_refine.Ubar, nurbs_refine.pu, 0, 1);
[Pedge, Sedge, Pmm, Smm, edgeMetaLeft] = assemble_vertical_edge_local( ...
    nurbs_refine.VBreaks, nurbs_refine.Vbar, nurbs_refine.pv, ...
    Uders_left(1, 1:2)', Uders_left(2, 1:2)', ...
    nurbs_refine.left_dofs_2_layers, -a, -1, ...
    pwx, pwy, alpha, Omega_area, a, inv_df, Ky);
P(nurbs_refine.left_dofs_2_layers, nurbs_refine.left_dofs_2_layers) = ...
    P(nurbs_refine.left_dofs_2_layers, nurbs_refine.left_dofs_2_layers) + Pedge.Ppp;
P(nurbs_refine.left_dofs_2_layers, pw_dofs) = ...
    P(nurbs_refine.left_dofs_2_layers, pw_dofs) + Pedge.Ppm;
P(pw_dofs, nurbs_refine.left_dofs_2_layers) = ...
    P(pw_dofs, nurbs_refine.left_dofs_2_layers) + Pedge.Pmp;
S(nurbs_refine.left_dofs_2_layers, nurbs_refine.left_dofs_2_layers) = ...
    S(nurbs_refine.left_dofs_2_layers, nurbs_refine.left_dofs_2_layers) + Sedge.Spp;
S(nurbs_refine.left_dofs_2_layers, pw_dofs) = ...
    S(nurbs_refine.left_dofs_2_layers, pw_dofs) + Sedge.Spm;
S(pw_dofs, nurbs_refine.left_dofs_2_layers) = ...
    S(pw_dofs, nurbs_refine.left_dofs_2_layers) + Sedge.Smp;
Ppw = Ppw + Pmm;
Spw = Spw + Smm;

% Right edge: u = 1, x = a, n = [1; 0]
Uders_right = bspbasisDers(nurbs_refine.Ubar, nurbs_refine.pu, 1, 1);
[Pedge, Sedge, Pmm, Smm, edgeMetaRight] = assemble_vertical_edge_local( ...
    nurbs_refine.VBreaks, nurbs_refine.Vbar, nurbs_refine.pv, ...
    Uders_right(1, end-1:end)', Uders_right(2, end-1:end)', ...
    nurbs_refine.right_dofs_2_layers, a, 1, ...
    pwx, pwy, alpha, Omega_area, a, inv_df, Ky);
P(nurbs_refine.right_dofs_2_layers, nurbs_refine.right_dofs_2_layers) = ...
    P(nurbs_refine.right_dofs_2_layers, nurbs_refine.right_dofs_2_layers) + Pedge.Ppp;
P(nurbs_refine.right_dofs_2_layers, pw_dofs) = ...
    P(nurbs_refine.right_dofs_2_layers, pw_dofs) + Pedge.Ppm;
P(pw_dofs, nurbs_refine.right_dofs_2_layers) = ...
    P(pw_dofs, nurbs_refine.right_dofs_2_layers) + Pedge.Pmp;
S(nurbs_refine.right_dofs_2_layers, nurbs_refine.right_dofs_2_layers) = ...
    S(nurbs_refine.right_dofs_2_layers, nurbs_refine.right_dofs_2_layers) + Sedge.Spp;
S(nurbs_refine.right_dofs_2_layers, pw_dofs) = ...
    S(nurbs_refine.right_dofs_2_layers, pw_dofs) + Sedge.Spm;
S(pw_dofs, nurbs_refine.right_dofs_2_layers) = ...
    S(pw_dofs, nurbs_refine.right_dofs_2_layers) + Sedge.Smp;
Ppw = Ppw + Pmm;
Spw = Spw + Smm;

meta.t_edges = toc(t_edges);

P(pw_dofs, pw_dofs) = P(pw_dofs, pw_dofs) + sparse(Ppw);
S(pw_dofs, pw_dofs) = S(pw_dofs, pw_dofs) + sparse(Spw);

meta.edge_bottom = edgeMetaBottom;
meta.edge_top = edgeMetaTop;
meta.edge_left = edgeMetaLeft;
meta.edge_right = edgeMetaRight;
meta.n_pw = n_pw;
meta.t_total = toc(t_total);

fprintf('[DG][square_fast] edge assembly total = %.4fs\n', meta.t_total);
end

function [Pedge, Sedge, Pmm, Smm, meta] = assemble_horizontal_edge_local( ...
breaks, knot, degree, fixed_vals, fixed_ders, plus_dofs, y_const, normal_y, ...
    pwx, pwy, alpha, Omega_area, a, inv_df, Kx)
% Assemble one horizontal DG interface edge.
% Prepare quadrature and sparse trace triplets.
[gp, gw] = grule(10 * degree + 5);
n_gp = numel(gp);
n_elem = numel(breaks) - 1;
nq = n_elem * n_gp;
n_pw = numel(pwx);
m = numel(plus_dofs) / 2;

nnz_est = nq * 2 * (degree + 1);
rows = zeros(nnz_est, 1);
cols = zeros(nnz_est, 1);
bvals = zeros(nnz_est, 1);
avals = zeros(nnz_est, 1);
wq = zeros(nq, 1);
xq = zeros(nq, 1);
normal_scale = normal_y * inv_df / 2;

cursor = 1;
qid = 1;

% Sample the IGA traces and normal derivatives.
for e = 1:n_elem
    a_param = breaks(e);
    b_param = breaks(e + 1);
    span = findspan(knot, degree, breaks(e));
    active = [span-degree:span, m + (span-degree:span)];

    param_q = ((b_param - a_param) * gp + a_param + b_param) / 2;
    jac_q = ((b_param - a_param) / 2) * gw * (2 * a);
    x_q = -a + 2 * a * param_q;

    for i = 1:n_gp
        ders = bspbasisDers(knot, degree, param_q(i), 1);
        nvar = ders(1, :).';

        basis_plus = (nvar * fixed_vals);
        basis_plus = basis_plus(:);

        normal_der = nvar * fixed_ders;
        avg_plus = (normal_der(:) * normal_scale).';

        local_n = numel(active);
        rows(cursor:cursor+local_n-1) = qid;
        cols(cursor:cursor+local_n-1) = active;
        bvals(cursor:cursor+local_n-1) = basis_plus;
        avals(cursor:cursor+local_n-1) = avg_plus;
        cursor = cursor + local_n;

        wq(qid) = jac_q(i);
        xq(qid) = x_q(i);
        qid = qid + 1;
    end
end

% Assemble the weighted IGA trace operators.
Bplus = sparse(rows, cols, bvals, nq, numel(plus_dofs));
Aplus = sparse(rows, cols, avals, nq, numel(plus_dofs));
W = spdiags(wq, 0, nq, nq);
WBplus = W * Bplus;
WAplus = W * Aplus;

% Evaluate the plane-wave traces and normal derivatives.
phase_y = exp(-1i * alpha * y_const * pwy(:)).';
Bminus = exp(-1i * alpha * (xq(:) * pwx(:).')) .* phase_y / sqrt(Omega_area);
mult = (1i * alpha * normal_y / 2) * pwy(:).';
WconjBminus = bsxfun(@times, conj(Bminus), wq);

% Contract the quadrature data into the four DG blocks.
Pedge = struct();
Pedge.Ppp = Bplus.' * WBplus;
Pedge.Ppm = -Bplus.' * WconjBminus;
Pedge.Pmp = -Bminus.' * WBplus;

Sedge = struct();
Sedge.Spp = Bplus.' * WAplus;
Sedge.Spm = Bplus.' * bsxfun(@times, WconjBminus, mult);
Sedge.Smp = -Bminus.' * WAplus;

phase_mat = (phase_y.' * conj(phase_y)) / Omega_area;
Pmm = phase_mat .* Kx;
Smm = -Pmm .* (ones(n_pw, 1) * mult);

meta = struct();
meta.n_elem = n_elem;
meta.n_gp = n_gp;
meta.nq = nq;
meta.plus_dofs = numel(plus_dofs);
end

function [Pedge, Sedge, Pmm, Smm, meta] = assemble_vertical_edge_local( ...
breaks, knot, degree, fixed_vals, fixed_ders, plus_dofs, x_const, normal_x, ...
    pwx, pwy, alpha, Omega_area, a, inv_df, Ky)
% Assemble one vertical DG interface edge.
% Prepare quadrature and sparse trace triplets.
[gp, gw] = grule(10 * degree + 5);
n_gp = numel(gp);
n_elem = numel(breaks) - 1;
nq = n_elem * n_gp;
n_pw = numel(pwx);

nnz_est = nq * 2 * (degree + 1);
rows = zeros(nnz_est, 1);
cols = zeros(nnz_est, 1);
bvals = zeros(nnz_est, 1);
avals = zeros(nnz_est, 1);
wq = zeros(nq, 1);
yq = zeros(nq, 1);
normal_scale = normal_x * inv_df / 2;

cursor = 1;
qid = 1;

% Sample the IGA traces and normal derivatives.
for e = 1:n_elem
    a_param = breaks(e);
    b_param = breaks(e + 1);
    span = findspan(knot, degree, breaks(e));
    active_var = span-degree:span;
    active = reshape([2 * active_var - 1; 2 * active_var], 1, []);

    param_q = ((b_param - a_param) * gp + a_param + b_param) / 2;
    jac_q = ((b_param - a_param) / 2) * gw * (2 * a);
    y_q = -a + 2 * a * param_q;

    for i = 1:n_gp
        ders = bspbasisDers(knot, degree, param_q(i), 1);
        nvar = ders(1, :);

        basis_plus = fixed_vals * nvar;
        basis_plus = basis_plus(:);

        normal_der = fixed_ders * nvar;
        avg_plus = (normal_der(:) * normal_scale).';

        local_n = numel(active);
        rows(cursor:cursor+local_n-1) = qid;
        cols(cursor:cursor+local_n-1) = active;
        bvals(cursor:cursor+local_n-1) = basis_plus;
        avals(cursor:cursor+local_n-1) = avg_plus;
        cursor = cursor + local_n;

        wq(qid) = jac_q(i);
        yq(qid) = y_q(i);
        qid = qid + 1;
    end
end

% Assemble the weighted IGA trace operators.
Bplus = sparse(rows, cols, bvals, nq, numel(plus_dofs));
Aplus = sparse(rows, cols, avals, nq, numel(plus_dofs));
W = spdiags(wq, 0, nq, nq);
WBplus = W * Bplus;
WAplus = W * Aplus;

% Evaluate the plane-wave traces and normal derivatives.
phase_x = exp(-1i * alpha * x_const * pwx(:)).';
Bminus = exp(-1i * alpha * (yq(:) * pwy(:).')) .* phase_x / sqrt(Omega_area);
mult = (1i * alpha * normal_x / 2) * pwx(:).';
WconjBminus = bsxfun(@times, conj(Bminus), wq);

% Contract the quadrature data into the four DG blocks.
Pedge = struct();
Pedge.Ppp = Bplus.' * WBplus;
Pedge.Ppm = -Bplus.' * WconjBminus;
Pedge.Pmp = -Bminus.' * WBplus;

Sedge = struct();
Sedge.Spp = Bplus.' * WAplus;
Sedge.Spm = Bplus.' * bsxfun(@times, WconjBminus, mult);
Sedge.Smp = -Bminus.' * WAplus;

phase_mat = (phase_x.' * conj(phase_x)) / Omega_area;
Pmm = phase_mat .* Ky;
Smm = -Pmm .* (ones(n_pw, 1) * mult);

meta = struct();
meta.n_elem = n_elem;
meta.n_gp = n_gp;
meta.nq = nq;
meta.plus_dofs = numel(plus_dofs);
end

function F = interval_ft_general_local(q, aL, aR, alpha)
% Integrate one Fourier mode over an arbitrary interval.
F = zeros(size(q));
idx0 = (q == 0);
F(idx0) = aR - aL;

qq = q(~idx0);
F(~idx0) = (exp(1i * alpha * qq * aR) - exp(1i * alpha * qq * aL)) ./ (1i * alpha * qq);
end
