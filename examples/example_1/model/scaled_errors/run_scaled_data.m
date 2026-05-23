function run_scaled_data()
%Run scaled data.

clc; close all; format short e;

activate_example_workflow('scaled_errors', ...
    {'nurbs', 'dg', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core', 'solver'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dataDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(dataDir);

% ---------------- user params ----------------
Example       = 'Example_1';
alpha_list    = [1.0, 1.5, 2.0, 2.5];
nElem_list    = [2 4 8 12 16];

% current runs: p = 1, 2
t_list        = [0, 1];     % p = 1+t -> [1,2]

n_eigenvalues = 1;

Nc_base       = 4;
beta          = 20;
n_gp          = 10;
inner_cheb_n  = 48;
pw_fft_grid_n = 256;

% ---------- existing reference file (fixed p=3) ----------
refRunMat      = fullfile(pwd, 'result', 'Example_1_verify49', ...
    'reference', 'p_3', 'refine_06_Nc_48', 'run.mat');
refNcFallback  = 48;
ref_pdeg       = 3;

% error evaluation grid
L              = 4;
a              = 0.2;
dx_in          = 1e-2;
dx_out         = 1e-2;
chunkSize      = 10000;

% PRIMME
primme_tol         = 1e-12;
primme_maxit       = 5e7;
primme_method      = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;
eps_diag           = 1e-12;
iface_reg          = 1e-12;

% save controls
save_eigenvectors = true;
save_nurbs        = true;
save_pw_index     = true;
use_pw_cache      = true;

% ---------------- result root ----------------
resultRoot = fullfile(pwd, 'result', 'scaled_errors');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

cacheRoot = fullfile(resultRoot, 'cache_pw');
if ~exist(cacheRoot, 'dir'), mkdir(cacheRoot); end

% ---------------- load existing reference ----------------
if ~exist(refRunMat, 'file')
    error('Reference file not found:\n%s', refRunMat);
end

fprintf('[LOAD REF] %s\n', refRunMat);
ref = load_run_data(refRunMat, refNcFallback);
ref.u  = ref.uh(:,1);
ref.uI = ref.u(1:ref.nNURBS);
ref.uA = ref.u(ref.nNURBS+1 : ref.nNURBS + size(ref.k_pw,1));

% ---------------- plan Nc from h ----------------
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

% loop over p = 1,2
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
        eL2      = [];
        eDG      = [];
        rhs49    = [];
        ratio49  = [];

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
                save_pw_index, use_pw_cache, cacheRoot, runDir, inner_cheb_n, pw_fft_grid_n);
            opts.mesh_mode = 'nElem';
            opts.nElem     = nElem;

            if ~exist(runMat, 'file')
                fprintf('\n------------------------------------------------------------\n');
                fprintf('[RUN ] p = %d, alpha = %.2f, nElem = %d, h = %.6f, Nc = %d\n', ...
                    pdeg, al, nElem, h, Nc);
                solve_iga_pw_dg(0, t, Nc, n_eigenvalues, opts);
            else
                fprintf('[LOAD] p = %d, alpha = %.2f, nElem = %d, h = %.6f, Nc = %d\n', ...
                    pdeg, al, nElem, h, Nc);
            end

            rr = load_run_data(runMat, Nc);

            u  = rr.uh(:,1);
            uI = u(1:rr.nNURBS);
            uA = u(rr.nNURBS+1 : rr.nNURBS + size(rr.k_pw,1));

            % phase alignment using common PW modes
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

            errL2 = compute_L2_error( ...
                rr.nurbs_refine, uI, uA, rr.k_pw, ...
                ref.nurbs_refine, ref.uI, ref.uA, ref.k_pw, ...
                L, a, dx_in, dx_out, chunkSize);

            errDG = compute_DG_error( ...
                rr.nurbs_refine, uI, uA, rr.k_pw, ...
                ref.nurbs_refine, ref.uI, ref.uA, ref.k_pw, ...
                L, a, dx_in, dx_out, chunkSize, sigma);

            rhs_now   = (rr.h ^ gamma49) * errDG;
            ratio_now = errL2 / max(rhs_now, eps);

            nElem_ok(end+1,1) = nElem; %#ok<AGROW>
            h_ok(end+1,1)     = rr.h;  %#ok<AGROW>
            Nc_ok(end+1,1)    = Nc;    %#ok<AGROW>
            dof_ok(end+1,1)   = rr.nDOF; %#ok<AGROW>
            lam_ok(end+1,1)   = real(rr.lambda(1)); %#ok<AGROW>
            eL2(end+1,1)      = errL2; %#ok<AGROW>
            eDG(end+1,1)      = errDG; %#ok<AGROW>
            rhs49(end+1,1)    = rhs_now; %#ok<AGROW>
            ratio49(end+1,1)  = ratio_now; %#ok<AGROW>

            fprintf('[ERR ] p=%d alpha=%.2f nElem=%d h=%.6f Nc=%d : L2=%.3e, DG=%.3e, h^g*DG=%.3e, ratio=%.3e\n', ...
                pdeg, al, nElem, rr.h, Nc, errL2, errDG, rhs_now, ratio_now);
        end

        T = table(repmat(al, numel(nElem_ok), 1), repmat(gamma49, numel(nElem_ok), 1), ...
            nElem_ok, h_ok, Nc_ok, dof_ok, lam_ok, eL2, eDG, rhs49, ratio49, ...
            'VariableNames', {'alpha','gamma49','nElem','h','Nc','dof','lambda1', ...
            'u1_L2','u1_DG','h_gamma_DG','ratio49'});

        csvFile = fullfile(caseRoot, 'errors.csv');
        matFile = fullfile(caseRoot, 'errors.mat');
        writetable(T, csvFile);
        save(matFile, 'T', 'al', 'gamma49', 'nElem_list', 'Nc_base', ...
            'refRunMat', 'refNcFallback', 'ref_pdeg', 'dx_in', 'dx_out', 'chunkSize');

        Tall = [Tall; T]; %#ok<AGROW>
    end

    summaryDir = fullfile(resultRoot, sprintf('p_%d', pdeg));
    if ~exist(summaryDir, 'dir'), mkdir(summaryDir); end
    allCsv = fullfile(summaryDir, 'summary.csv');
    writetable(Tall, allCsv);
end

end

function opts = make_opts(Example, beta, n_gp, primme_tol, primme_maxit, ...
    primme_method, primme_reportLevel, eps_diag, iface_reg, save_eigenvectors, save_nurbs, ...
    save_pw_index, use_pw_cache, cacheRoot, outDir, inner_cheb_n, pw_fft_grid_n)
%Build solver options.

opts = struct();
opts.Example            = Example;
opts.beta               = beta;
opts.n_gp               = n_gp;
opts.inner_cheb_n       = inner_cheb_n;
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

opts.save_matrices      = false;
opts.save_mat           = false;
end

function Nc = choose_Nc_from_h(h, h0, Nc0, alpha)
%Choose the plane-wave cutoff from the mesh size.
Nc = ceil(Nc0 * (h0 / h)^(1 / alpha));
Nc = max(2, Nc);
end

function tag = alpha_to_tag(alpha)
%Compute to tag.
tag = strrep(sprintf('%.2f', alpha), '.', 'p');
end

function rr = load_run_data(runMat, NcFallback)
%Load run data.
S = load(runMat);
if isfield(S, 'run')
    R = S.run;
else
    R = S;
end

if ~isfield(R, 'uh') || isempty(R.uh)
    error('run.mat does not contain eigenvectors uh:\n%s', runMat);
end
rr.uh = R.uh;

if ~isfield(R, 'lambda') || isempty(R.lambda)
    error('run.mat does not contain lambda:\n%s', runMat);
end
rr.lambda = R.lambda(:).';

if isfield(R, 'k_pw') && ~isempty(R.k_pw)
    rr.k_pw = R.k_pw;
else
    [rr.k_pw, ~] = generate_p_vec(NcFallback);
end

if isfield(R, 'n_dofs_total') && ~isempty(R.n_dofs_total)
    rr.nDOF = double(R.n_dofs_total);
elseif isfield(R, 'meta') && isfield(R.meta, 'n_dofs_total')
    rr.nDOF = double(R.meta.n_dofs_total);
else
    rr.nDOF = size(R.uh, 1);
end

if isfield(R, 'n_dofs_nurbs') && ~isempty(R.n_dofs_nurbs)
    rr.nNURBS = double(R.n_dofs_nurbs);
elseif isfield(R, 'meta') && isfield(R.meta, 'n_dofs_nurbs')
    rr.nNURBS = double(R.meta.n_dofs_nurbs);
else
    rr.nNURBS = size(R.uh,1) - size(rr.k_pw,1);
end

if ~isfield(R, 'nurbs_refine') || isempty(R.nurbs_refine)
    error('run.mat does not contain nurbs_refine:\n%s', runMat);
end
rr.nurbs_refine = R.nurbs_refine;

if isfield(R, 'nurbs_original')
    rr.nurbs_original = R.nurbs_original;
else
    rr.nurbs_original = [];
end

if isfield(R, 'meta') && isfield(R.meta, 'h') && ~isempty(R.meta.h)
    rr.h = R.meta.h;
else
    rr.h = 2 * 0.2 * estimate_h_parametric(rr.nurbs_refine);
end
end

function errL2 = compute_L2_error( ...
nurbs_curr, uI_curr, uA_curr, k_curr, ...
    nurbs_ref,  uI_ref,  uA_ref,  k_ref, ...
    L, a, dx_in, dx_out, chunkSize)

% Compute the L2 error on the comparison grid.
[xi, yi, wA_in] = grid_points_square(-a, a, dx_in);
vI_curr = iga_eval_val(nurbs_curr, uI_curr, xi, yi, a);
vI_ref  = iga_eval_val(nurbs_ref,  uI_ref,  xi, yi, a);
err_in  = sum(abs(vI_curr - vI_ref).^2) * wA_in;

[xo, yo, wA_out] = grid_points_outer(L, a, dx_out);
err_out = 0;
nPts = numel(xo);
k = 1;
while k <= nPts
    k2 = min(nPts, k + chunkSize - 1);
    X  = xo(k:k2);
    Y  = yo(k:k2);

    vA_curr = pw_eval_val(uA_curr, k_curr, X, Y, L);
    vA_ref  = pw_eval_val(uA_ref,  k_ref,  X, Y, L);

    err_out = err_out + sum(abs(vA_curr - vA_ref).^2);
    k = k2 + 1;
end
err_out = err_out * wA_out;

errL2 = sqrt(err_in + err_out);
end

function errDG = compute_DG_error( ...
nurbs_curr, uI_curr, uA_curr, k_curr, ...
    nurbs_ref,  uI_ref,  uA_ref,  k_ref, ...
    L, a, dx_in, dx_out, chunkSize, sigma)

% Compute the DG error on the comparison grid.
[xi, yi, wA_in] = grid_points_square(-a, a, dx_in);
[vI_curr, gxI_curr, gyI_curr] = iga_eval_val_grad(nurbs_curr, uI_curr, xi, yi, a);
[vI_ref,  gxI_ref,  gyI_ref ] = iga_eval_val_grad(nurbs_ref,  uI_ref,  xi, yi, a);

de  = vI_curr - vI_ref;
dex = gxI_curr - gxI_ref;
dey = gyI_curr - gyI_ref;
H1_in = sum(abs(de).^2 + abs(dex).^2 + abs(dey).^2) * wA_in;

[xo, yo, wA_out] = grid_points_outer(L, a, dx_out);
H1_out = 0;
nPts = numel(xo);
k = 1;
while k <= nPts
    k2 = min(nPts, k + chunkSize - 1);
    X  = xo(k:k2);
    Y  = yo(k:k2);

    [vA_curr, gxA_curr, gyA_curr] = pw_eval_val_grad(uA_curr, k_curr, X, Y, L);
    [vA_ref,  gxA_ref,  gyA_ref ] = pw_eval_val_grad(uA_ref,  k_ref,  X, Y, L);

    de  = vA_curr - vA_ref;
    dex = gxA_curr - gxA_ref;
    dey = gyA_curr - gyA_ref;
    H1_out = H1_out + sum(abs(de).^2 + abs(dex).^2 + abs(dey).^2);
    k = k2 + 1;
end
H1_out = H1_out * wA_out;

[xg, yg, wL] = boundary_points_square(a, dx_in);
vA_curr = pw_eval_val(uA_curr, k_curr, xg, yg, L);
vA_ref  = pw_eval_val(uA_ref,  k_ref,  xg, yg, L);
vI_curr = iga_eval_val(nurbs_curr, uI_curr, xg, yg, a);
vI_ref  = iga_eval_val(nurbs_ref,  uI_ref,  xg, yg, a);

jump = (vA_curr - vA_ref) - (vI_curr - vI_ref);
J2 = sum(abs(jump).^2) * wL;

errDG = sqrt(H1_in + H1_out + sigma * J2);
end

function val = pw_eval_val(coeff, p_vec, X, Y, L)
%Evaluate the field value.
F = [X(:)'; Y(:)'];
expo = exp((1i * 2*pi / L) * (p_vec * F));
val = (coeff.' * expo) / L;
val = val(:);
end

function [val, gx, gy] = pw_eval_val_grad(coeff, p_vec, X, Y, L)
%Evaluate the field value and gradient.
F = [X(:)'; Y(:)'];
expo = exp((1i * 2*pi / L) * (p_vec * F));

val = (coeff.' * expo) / L;
fac = (1i * 2*pi / L) / L;
gx  = ((coeff .* p_vec(:,1)).' * expo) * fac;
gy  = ((coeff .* p_vec(:,2)).' * expo) * fac;

val = val(:);
gx  = gx(:);
gy  = gy(:);
end

function val = iga_eval_val(nurbs, coeff, X, Y, a)
%Evaluate the field value.
[val, ~, ~] = iga_eval_val_grad(nurbs, coeff, X, Y, a);
end

function [val, gx, gy] = iga_eval_val_grad(nurbs, coeff, X, Y, a)
%Evaluate the field value and gradient.
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
%Evaluate basis values and first derivatives.
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
%Locate an index or object used by the computation.
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
%Build midpoint quadrature points in the square.
x = xmin + dx/2 : dx : xmax - dx/2;
[Xg, Yg] = meshgrid(x, x);
X = Xg(:);
Y = Yg(:);
wA = dx * dx;
end

function [X, Y, wA] = grid_points_outer(L, a, dx)
%Build midpoint quadrature points outside the square.
x = -L/2 + dx/2 : dx : L/2 - dx/2;
[Xg, Yg] = meshgrid(x, x);
mask = ~(Xg >= -a & Xg <= a & Yg >= -a & Yg <= a);
X = Xg(mask);
Y = Yg(mask);
wA = dx * dx;
end

function [X, Y, wL] = boundary_points_square(a, ds)
%Build midpoint quadrature points on the square boundary.
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
%Estimate the parametric mesh size.
Uu = unique(nurbs_refine.Ubar(:).');
Vv = unique(nurbs_refine.Vbar(:).');
h = max(max(diff(Uu)), max(diff(Vv)));
end

function [p_vec, n_pw_basis] = generate_p_vec(Nc)
%Generate p vec.
N = floor(Nc);
p_vec = zeros((2*N+1)^2, 2);
n_pw_basis = 0;

for ii = -N:N
    m = floor(sqrt(N^2 - ii^2));
    for jj = -m:m
        n_pw_basis = n_pw_basis + 1;
        p_vec(n_pw_basis,:) = [ii, jj];
    end
end

p_vec = p_vec(1:n_pw_basis, :);
end
