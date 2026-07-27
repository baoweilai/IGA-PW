function run_scaled_data()
% Generate Example 1 scaled-error data.

clc; close all; format short e;

% Set workflow paths and the data directory.
activate_example_workflow('scaled_errors', ...
    {'nurbs', 'dg', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dataDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(dataDir);

% Set the scaling and discretization parameters.
Example       = 'Example_1';
alpha_list    = [1.0, 2.0];
nElem_list    = [2 4 8 12 16];

% Use degree p = 2.
t_list        = 1;

n_eigenvalues = 4;

Nc_base       = 4;
beta          = 20;
n_gp          = 10;
inner_cheb_n  = 100;
inner_quad_n  = 400;
pw_fft_grid_n = 500;

% Set the saved four-vector reference.
refRunMat      = fullfile(pwd, 'result', 'Example_1', 'Nc_45', ...
    'p_3', 'refine_08', 'run.mat');
ref_pdeg       = 3;

% Set the error quadrature grid.
L              = 4;
a              = 0.2;
dx_in          = 1e-2;
dx_out         = 1e-2;
chunkSize      = 10000;

% Set the eigensolver parameters.
primme_tol         = 1e-12;
primme_maxit       = 5e7;
primme_method      = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;
eps_diag           = 1e-12;
iface_reg          = 1e-12;

% Set the saved run fields.
save_eigenvectors = true;
save_nurbs        = true;
save_pw_index     = true;
use_pw_cache      = true;

% Create the result and cache directories.
resultRoot = fullfile(pwd, 'result', 'scaled_errors');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

cacheRoot = fullfile(resultRoot, 'cache_pw');
if ~exist(cacheRoot, 'dir'), mkdir(cacheRoot); end

% Load the reference eigenvectors.
if ~exist(refRunMat, 'file')
    error('Reference file not found:\n%s', refRunMat);
end

fprintf('[LOAD REF] %s\n', refRunMat);
ref = load_run_data(refRunMat);
assert(size(ref.uh, 2) >= n_eigenvalues, ...
    'Reference file must contain the first four eigenvectors.');
ref.u  = ref.uh(:,1);
ref.uI = ref.u(1:ref.nNURBS);
ref.uA = ref.u(ref.nNURBS+1 : ref.nNURBS + size(ref.k_pw,1));

% Compute the cutoff associated with each mesh size.
h0 = 2 * a / nElem_list(1);

fprintf('================ plan (nElem, h, Nc) ================\n');
for ia = 1:numel(alpha_list)
    al = alpha_list(ia);
    fprintf('alpha = %.2f\n', al);
    for ir = 1:numel(nElem_list)
        nElem = nElem_list(ir);
        h = 2 * a / nElem;
        Nc = choose_Nc_from_h(h, h0, Nc_base, al);
        fprintf('  nElem = %d, h = %.6f, Nc = %d\n', nElem, h, Nc);
    end
end

% Solve the configured degree and scaling cases.
for it = 1:numel(t_list)
    t = t_list(it);
    pdeg = 1 + t;

    fprintf('\n############################################################\n');
    fprintf('[RUN BLOCK] current p = %d, reference p = %d\n', pdeg, ref_pdeg);
    fprintf('############################################################\n');

    Tall = table();

    for ia = 1:numel(alpha_list)
        al = alpha_list(ia);
        gamma49 = (3 - al) / (2 * al);
        alphaTag = alpha_to_tag(al);

        caseRoot = fullfile(resultRoot, sprintf('alpha_%s', alphaTag), sprintf('p_%d', pdeg));
        if ~exist(caseRoot, 'dir'), mkdir(caseRoot); end

        nElem_ok = [];
        h_ok     = [];
        Nc_ok    = [];
        dof_ok   = [];
        lam_ok   = [];
        lam3_ok  = [];
        lam4_ok  = [];
        eL2      = [];
        eDG      = [];
        eL2_u3   = [];
        eDG_u3   = [];
        rhs49    = [];
        ratio49  = [];
        rhs49_u3 = [];
        ratio49_u3 = [];

        fprintf('\n============================================================\n');
        fprintf('[ALPHA] alpha = %.2f, gamma = %.6f, current p = %d\n', al, gamma49, pdeg);

        for ir = 1:numel(nElem_list)
            nElem = nElem_list(ir);
            h = 2 * a / nElem;
            Nc = choose_Nc_from_h(h, h0, Nc_base, al);

            runDir = fullfile(caseRoot, sprintf('nElem_%03d_Nc_%02d', nElem, Nc));
            if ~exist(runDir, 'dir'), mkdir(runDir); end
            runMat = fullfile(runDir, 'run.mat');

            opts = make_opts(Example, beta, n_gp, primme_tol, primme_maxit, ...
                primme_method, primme_reportLevel, eps_diag, iface_reg, save_eigenvectors, save_nurbs, ...
                save_pw_index, use_pw_cache, cacheRoot, runDir, inner_cheb_n, inner_quad_n, pw_fft_grid_n);
            opts.mesh_mode = 'nElem';
            opts.nElem     = nElem;

            if ~run_matches_parameters(runMat, n_eigenvalues, inner_cheb_n, ...
                    inner_quad_n, pw_fft_grid_n, pdeg, nElem, Nc)
                fprintf('\n------------------------------------------------------------\n');
                fprintf('[RUN ] p = %d, alpha = %.2f, nElem = %d, h = %.6f, Nc = %d\n', ...
                    pdeg, al, nElem, h, Nc);
                solve_iga_pw_dg(0, t, Nc, n_eigenvalues, opts);
            else
                fprintf('[LOAD] p = %d, alpha = %.2f, nElem = %d, h = %.6f, Nc = %d\n', ...
                    pdeg, al, nElem, h, Nc);
            end

            rr = load_run_data(runMat);
            assert(size(rr.uh, 2) >= n_eigenvalues, ...
                'Run file contains fewer than four eigenvectors: %s', runMat);

            u  = rr.uh(:,1);
            uI = u(1:rr.nNURBS);
            uA = u(rr.nNURBS+1 : rr.nNURBS + size(rr.k_pw,1));

            % Align the first eigenvector using common plane-wave modes.
            [tf, loc] = ismember(rr.k_pw, ref.k_pw, 'rows');
            if ~all(tf)
                error('Current PW modes are not contained in reference PW modes.');
            end

            uA_pad = zeros(size(ref.uA));
            uA_pad(loc(tf)) = uA(tf);

            alpha_phase = ref.uA' * uA_pad;
            if abs(alpha_phase) > 1e-14
                phase = exp(-1i * angle(alpha_phase));
            else
                phase = 1;
            end

            uI = uI * phase;
            uA = uA * phase;

            sigma = beta * (Nc + 1 / max(rr.h, 1e-14));

            [uI34, uA34, refI34, refA34] = align_eigenspace34( ...
                rr, ref, L, a, dx_in, dx_out, chunkSize);

            [errL2_pair, errDG_pair] = compute_vector_errors( ...
                rr.nurbs_refine, [uI, uI34(:,1)], ...
                [uA, uA34(:,1)], rr.k_pw, ...
                ref.nurbs_refine, [ref.uI, refI34(:,1)], ...
                [ref.uA, refA34(:,1)], ref.k_pw, ...
                L, a, dx_in, dx_out, chunkSize, sigma);

            errL2 = errL2_pair(1);
            errDG = errDG_pair(1);
            errL2_u3_now = errL2_pair(2);
            errDG_u3_now = errDG_pair(2);

            rhs_now   = (rr.h ^ gamma49) * errDG;
            ratio_now = errL2 / max(rhs_now, eps);
            rhs_u3_now = (rr.h ^ gamma49) * errDG_u3_now;
            ratio_u3_now = errL2_u3_now / max(rhs_u3_now, eps);

            nElem_ok(end+1,1) = nElem; %#ok<AGROW>
            h_ok(end+1,1)     = rr.h;  %#ok<AGROW>
            Nc_ok(end+1,1)    = Nc;    %#ok<AGROW>
            dof_ok(end+1,1)   = rr.nDOF; %#ok<AGROW>
            lam_ok(end+1,1)   = real(rr.lambda(1)); %#ok<AGROW>
            lam3_ok(end+1,1)  = real(rr.lambda(3)); %#ok<AGROW>
            lam4_ok(end+1,1)  = real(rr.lambda(4)); %#ok<AGROW>
            eL2(end+1,1)      = errL2; %#ok<AGROW>
            eDG(end+1,1)      = errDG; %#ok<AGROW>
            eL2_u3(end+1,1)   = errL2_u3_now; %#ok<AGROW>
            eDG_u3(end+1,1)   = errDG_u3_now; %#ok<AGROW>
            rhs49(end+1,1)    = rhs_now; %#ok<AGROW>
            ratio49(end+1,1)  = ratio_now; %#ok<AGROW>
            rhs49_u3(end+1,1) = rhs_u3_now; %#ok<AGROW>
            ratio49_u3(end+1,1) = ratio_u3_now; %#ok<AGROW>

            fprintf(['[ERR ] p=%d alpha=%.2f nElem=%d h=%.6f Nc=%d : ' ...
                'u1[L2=%.3e DG=%.3e h^gDG=%.3e] ' ...
                'u3[L2=%.3e DG=%.3e h^gDG=%.3e]\n'], ...
                pdeg, al, nElem, rr.h, Nc, errL2, errDG, rhs_now, ...
                errL2_u3_now, errDG_u3_now, rhs_u3_now);
        end

        % Save the errors for this scaling exponent.
        T = table(repmat(al, numel(nElem_ok), 1), repmat(gamma49, numel(nElem_ok), 1), ...
            nElem_ok, h_ok, Nc_ok, dof_ok, lam_ok, lam3_ok, lam4_ok, ...
            eL2, eDG, rhs49, ratio49, eL2_u3, eDG_u3, rhs49_u3, ratio49_u3, ...
            'VariableNames', {'alpha','gamma49','nElem','h','Nc','dof', ...
            'lambda1','lambda3','lambda4','u1_L2','u1_DG','h_gamma_DG', ...
            'ratio49','u3_L2','u3_DG','u3_h_gamma_DG','u3_ratio49'});

        csvFile = fullfile(caseRoot, 'errors.csv');
        matFile = fullfile(caseRoot, 'errors.mat');
        writetable(T, csvFile);
        save(matFile, 'T', 'al', 'gamma49', 'nElem_list', 'Nc_base', ...
            'refRunMat', 'ref_pdeg', 'dx_in', 'dx_out', ...
            'chunkSize', 'n_eigenvalues', 'inner_cheb_n', 'inner_quad_n', ...
            'pw_fft_grid_n');

        Tall = [Tall; T]; %#ok<AGROW>
    end

    % Save the combined degree summary.
    summaryDir = fullfile(resultRoot, sprintf('p_%d', pdeg));
    if ~exist(summaryDir, 'dir'), mkdir(summaryDir); end
    allCsv = fullfile(summaryDir, 'summary.csv');
    writetable(Tall, allCsv);
end

end

function opts = make_opts(Example, beta, n_gp, primme_tol, primme_maxit, ...
    primme_method, primme_reportLevel, eps_diag, iface_reg, save_eigenvectors, save_nurbs, ...
    save_pw_index, use_pw_cache, cacheRoot, outDir, inner_cheb_n, inner_quad_n, pw_fft_grid_n)
% Build solver options.

opts = struct();
opts.Example            = Example;
opts.beta               = beta;
opts.n_gp               = n_gp;
opts.inner_cheb_n       = inner_cheb_n;
opts.inner_quad_n       = inner_quad_n;
opts.pw_fft_grid_n      = pw_fft_grid_n;

opts.primme_tol         = primme_tol;
opts.primme_maxit       = primme_maxit;
opts.primme_method      = primme_method;
opts.primme_reportLevel = primme_reportLevel;
opts.eps_diag           = eps_diag;
opts.iface_reg          = iface_reg;

opts.save_eigenvectors  = save_eigenvectors;
opts.save_nurbs         = save_nurbs;
opts.save_pw_index      = save_pw_index;

opts.use_pw_cache       = use_pw_cache;
opts.cacheRoot          = cacheRoot;
opts.outDir             = outDir;
end

function Nc = choose_Nc_from_h(h, h0, Nc0, alpha)
% Choose the plane-wave cutoff from the mesh size.
Nc = ceil(Nc0 * (h0 / h)^(1 / alpha));
Nc = max(2, Nc);
end

function tag = alpha_to_tag(alpha)
% Encode alpha in a filename-safe tag.
tag = strrep(sprintf('%.2f', alpha), '.', 'p');
end

function rr = load_run_data(runMat)
% Load one saved scaled-error run.
S = load(runMat);
assert(isfield(S, 'run'), 'MAT file does not contain run: %s', runMat);
R = S.run;
assert(isfield(R, 'uh') && ~isempty(R.uh), ...
    'run.mat does not contain eigenvectors uh: %s', runMat);
rr.uh = R.uh;
assert(isfield(R, 'lambda') && ~isempty(R.lambda), ...
    'run.mat does not contain lambda: %s', runMat);
rr.lambda = R.lambda(:).';
assert(isfield(R, 'k_pw') && ~isempty(R.k_pw), ...
    'run.mat does not contain k_pw: %s', runMat);
rr.k_pw = R.k_pw;
assert(isfield(R, 'n_dofs_total') && ~isempty(R.n_dofs_total), ...
    'run.mat does not contain n_dofs_total: %s', runMat);
rr.nDOF = double(R.n_dofs_total);
assert(isfield(R, 'n_dofs_nurbs') && ~isempty(R.n_dofs_nurbs), ...
    'run.mat does not contain n_dofs_nurbs: %s', runMat);
rr.nNURBS = double(R.n_dofs_nurbs);
assert(isfield(R, 'nurbs_refine') && ~isempty(R.nurbs_refine), ...
    'run.mat does not contain nurbs_refine: %s', runMat);
rr.nurbs_refine = R.nurbs_refine;
assert(isfield(R, 'nurbs_original') && ~isempty(R.nurbs_original), ...
    'run.mat does not contain nurbs_original: %s', runMat);
rr.nurbs_original = R.nurbs_original;

if isfield(R, 'meta') && isfield(R.meta, 'h') && ~isempty(R.meta.h)
    rr.h = R.meta.h;
else
    rr.h = 2 * 0.2 * estimate_h_parametric(rr.nurbs_refine);
end
end

function tf = run_matches_parameters(runMat, nEigenvalues, innerChebN, ...
    innerQuadN, fftGridN, pdeg, nElem, Nc)
% Accept only data generated with the requested eigensolver and PW settings.
tf = false;
if exist(runMat, 'file') ~= 2
    return;
end

S = load(runMat, 'run');
if ~isfield(S, 'run') || ~isfield(S.run, 'meta') || ...
        ~isfield(S.run, 'uh') || size(S.run.uh, 2) < nEigenvalues
    return;
end

meta = S.run.meta;
required = {'inner_cheb_n', 'inner_quad_n', 'pw_fft_grid_n', ...
    'pu', 'nElem', 'Nc'};
if ~all(isfield(meta, required))
    return;
end

tf = meta.inner_cheb_n == innerChebN ...
    && meta.inner_quad_n == innerQuadN ...
    && meta.pw_fft_grid_n == fftGridN ...
    && meta.pu == pdeg ...
    && meta.nElem == nElem ...
    && meta.Nc == Nc;
end

function [uI34, uA34, refI34, refA34] = align_eigenspace34( ...
    curr, ref, L, a, dx_in, dx_out, chunkSize)
% L2-Procrustes alignment of the repeated third/fourth eigenspaces.
currI = curr.uh(1:curr.nNURBS, 3:4);
currA = curr.uh(curr.nNURBS + 1:curr.nNURBS + size(curr.k_pw, 1), 3:4);
refI = ref.uh(1:ref.nNURBS, 3:4);
refA = ref.uh(ref.nNURBS + 1:ref.nNURBS + size(ref.k_pw, 1), 3:4);

[xi, yi, wA_in] = grid_points_square(-a, a, dx_in);
vCurr = iga_eval_block(curr.nurbs_refine, currI, xi, yi, a);
vRef = iga_eval_block(ref.nurbs_refine, refI, xi, yi, a);
Mrr = vRef' * vRef * wA_in;
Mhh = vCurr' * vCurr * wA_in;
Mrh = vRef' * vCurr * wA_in;

[xo, yo, wA_out] = grid_points_outer(L, a, dx_out);
for k = 1:chunkSize:numel(xo)
    idx = k:min(numel(xo), k + chunkSize - 1);
    vCurr = pw_eval_block(currA, curr.k_pw, xo(idx), yo(idx), L);
    vRef = pw_eval_block(refA, ref.k_pw, xo(idx), yo(idx), L);
    Mrr = Mrr + vRef' * vRef * wA_out;
    Mhh = Mhh + vCurr' * vCurr * wA_out;
    Mrh = Mrh + vRef' * vCurr * wA_out;
end

Rr = chol(hermitian(Mrr));
Rh = chol(hermitian(Mhh));
C = Rr' \ (Mrh / Rh);
[left, ~, right] = svd(C, 'econ');
Q = right * left';

refI34 = refI / Rr;
refA34 = refA / Rr;
uI34 = (currI / Rh) * Q;
uA34 = (currA / Rh) * Q;
end

function [errL2, errDG] = compute_vector_errors( ...
    nurbs_curr, uI_curr, uA_curr, k_curr, ...
    nurbs_ref, uI_ref, uA_ref, k_ref, ...
    L, a, dx_in, dx_out, chunkSize, sigma)
% Compute L2 and broken-H1-plus-jump errors for several aligned vectors.
[xi, yi, wA_in] = grid_points_square(-a, a, dx_in);
[vCurr, gxCurr, gyCurr] = iga_eval_block(nurbs_curr, uI_curr, xi, yi, a);
[vRef, gxRef, gyRef] = iga_eval_block(nurbs_ref, uI_ref, xi, yi, a);

dv = vCurr - vRef;
dx = gxCurr - gxRef;
dy = gyCurr - gyRef;
L2sq = sum(abs(dv).^2, 1) * wA_in;
DGsq = sum(abs(dv).^2 + abs(dx).^2 + abs(dy).^2, 1) * wA_in;

[xo, yo, wA_out] = grid_points_outer(L, a, dx_out);
for k = 1:chunkSize:numel(xo)
    idx = k:min(numel(xo), k + chunkSize - 1);
    [vCurr, gxCurr, gyCurr] = pw_eval_block( ...
        uA_curr, k_curr, xo(idx), yo(idx), L);
    [vRef, gxRef, gyRef] = pw_eval_block( ...
        uA_ref, k_ref, xo(idx), yo(idx), L);
    dv = vCurr - vRef;
    dx = gxCurr - gxRef;
    dy = gyCurr - gyRef;
    L2sq = L2sq + sum(abs(dv).^2, 1) * wA_out;
    DGsq = DGsq + sum(abs(dv).^2 + abs(dx).^2 + abs(dy).^2, 1) * wA_out;
end

[xg, yg, wL] = boundary_points_square(a, dx_in);
vACurr = pw_eval_block(uA_curr, k_curr, xg, yg, L);
vARef = pw_eval_block(uA_ref, k_ref, xg, yg, L);
vICurr = iga_eval_block(nurbs_curr, uI_curr, xg, yg, a);
vIRef = iga_eval_block(nurbs_ref, uI_ref, xg, yg, a);
jump = (vACurr - vARef) - (vICurr - vIRef);
DGsq = DGsq + sigma * sum(abs(jump).^2, 1) * wL;

errL2 = sqrt(max(real(L2sq), 0));
errDG = sqrt(max(real(DGsq), 0));
end

function [val, gx, gy] = pw_eval_block(coeff, p_vec, X, Y, L)
% Evaluate one or more plane-wave fields and their gradients.
F = [X(:)'; Y(:)'];
expo = exp((1i * 2*pi / L) * (p_vec * F));
fac = (1i * 2*pi / L) / L;
val = ((coeff.' * expo) / L).';
gx = (((coeff .* p_vec(:,1)).' * expo) * fac).';
gy = (((coeff .* p_vec(:,2)).' * expo) * fac).';
end

function [val, gx, gy] = iga_eval_block(nurbs, coeff, X, Y, a)
% Evaluate one or more IGA fields and their gradients.
nFields = size(coeff, 2);
val = zeros(numel(X), nFields);
gx = zeros(numel(X), nFields);
gy = zeros(numel(X), nFields);
for j = 1:nFields
    [val(:,j), gx(:,j), gy(:,j)] = iga_eval_val_grad( ...
        nurbs, coeff(:,j), X, Y, a);
end
end

function A = hermitian(A)
% Remove quadrature-level loss of Hermitian symmetry.
A = (A + A') / 2;
end

function [val, gx, gy] = iga_eval_val_grad(nurbs, coeff, X, Y, a)
% Evaluate the field value and gradient.
pu = nurbs.pu;
pv = nurbs.pv;
U  = nurbs.Ubar(:).';
V  = nurbs.Vbar(:).';

mU = length(U) - pu - 1;
nV = length(V) - pv - 1;

val = zeros(numel(X),1);
gx  = zeros(numel(X),1);
gy  = zeros(numel(X),1);

for k = 1:numel(X)
    u = (X(k) + a) / (2*a);
    v = (Y(k) + a) / (2*a);
    u = max(0, min(1, u));
    v = max(0, min(1, v));

    spanU = findspan_local_local(mU-1, pu, u, U);
    spanV = findspan_local_local(nV-1, pv, v, V);

    [Nu, dNu] = bspline_basis_and_der1_local(U, pu, u, spanU);
    [Nv, dNv] = bspline_basis_and_der1_local(V, pv, v, spanV);

    s  = 0;
    su = 0;
    sv = 0;

    for j = (spanV-pv):spanV
        lv = j - (spanV-pv) + 1;
        for i = (spanU-pu):spanU
            lu = i - (spanU-pu) + 1;
            row = i + (j-1) * mU;
            c = coeff(row);

            s  = s  + c * Nu(lu)  * Nv(lv);
            su = su + c * dNu(lu) * Nv(lv);
            sv = sv + c * Nu(lu)  * dNv(lv);
        end
    end

    val(k) = s;
    gx(k)  = su / (2*a);
    gy(k)  = sv / (2*a);
end
end

function [N, dN] = bspline_basis_and_der1_local(U, p, u, span)
% Evaluate basis values and first derivatives.
ndu = zeros(p+1, p+1);
left = zeros(1, p+1);
right = zeros(1, p+1);

ndu(1,1) = 1.0;
for j = 1:p
    left(j+1)  = u - U(span+1-j);
    right(j+1) = U(span+j) - u;
    saved = 0.0;
    for r = 0:j-1
        ndu(j+1, r+1) = right(r+2) + left(j-r+1);
        temp = ndu(r+1, j) / ndu(j+1, r+1);
        ndu(r+1, j+1) = saved + right(r+2) * temp;
        saved = left(j-r+1) * temp;
    end
    ndu(j+1, j+1) = saved;
end

N = ndu(1:p+1, p+1).';

ders1 = zeros(1, p+1);
for r = 0:p
    d = 0.0;
    pk = p - 1;
    if r >= 1
        d = d + ndu(r, pk+1) / ndu(pk+2, r);
    end
    if r <= p-1
        d = d - ndu(r+1, pk+1) / ndu(pk+2, r+1);
    end
    ders1(r+1) = d;
end
dN = ders1 * p;
end

function span = findspan_local_local(n, p, u, U)
% Find the active knot span for a parameter value.
if u >= U(n+2)
    span = n + 1;
    return;
end
if u <= U(p+1)
    span = p + 1;
    return;
end

low = p + 1;
high = n + 2;
mid = floor((low + high) / 2);

while (u < U(mid) || u >= U(mid+1))
    if u < U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high) / 2);
end

span = mid;
end

function [X, Y, wA] = grid_points_square(xmin, xmax, dx)
% Build midpoint quadrature points in the square.
x = xmin + dx/2 : dx : xmax - dx/2;
[Xg, Yg] = meshgrid(x, x);
X = Xg(:);
Y = Yg(:);
wA = dx * dx;
end

function [X, Y, wA] = grid_points_outer(L, a, dx)
% Build midpoint quadrature points outside the square.
x = -L/2 + dx/2 : dx : L/2 - dx/2;
[Xg, Yg] = meshgrid(x, x);
mask = ~(Xg >= -a & Xg <= a & Yg >= -a & Yg <= a);
X = Xg(mask);
Y = Yg(mask);
wA = dx * dx;
end

function [X, Y, wL] = boundary_points_square(a, ds)
% Build midpoint quadrature points on the square boundary.
t = -a + ds/2 : ds : a - ds/2;

xb = t;  yb = -a*ones(size(t));
xt = t;  yt =  a*ones(size(t));
xl = -a*ones(size(t)); yl = t;
xr =  a*ones(size(t)); yr = t;

X = [xb, xt, xl, xr].';
Y = [yb, yt, yl, yr].';
wL = ds;
end

function h = estimate_h_parametric(nurbs_refine)
% Estimate the parametric mesh size.
Uu = unique(nurbs_refine.Ubar(:).');
Vv = unique(nurbs_refine.Vbar(:).');
h = max(max(diff(Uu)), max(diff(Vv)));
end
