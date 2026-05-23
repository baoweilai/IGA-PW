function summary = plot_preconditioner(runCfg)
%Plot preconditioner comparison data.

clc; close all;

assert(exist('runCfg', 'var') == 1, 'plot_preconditioner requires runCfg.');

scriptDir = fileparts(mfilename('fullpath'));
exampleDir = fileparts(scriptDir);
activate_example_workflow('preconditioner', ...
    {'nurbs', 'iga', 'assembly', 'operators', 'core', 'solver'});

defaultRunCfg = struct( ...
    'Nc', 20, ...
    'p', 1, ...
    'refine_list', 1:7, ...
    'n_eigenvalues', 1, ...
    'beta', 20, ...
    'n_gp', 10, ...
    'primme_tol', 1e-8, ...
    'primme_reportLevel', 2, ...
    'scf_tol_lambda', 1e-6, ...
    'scf_maxit', 60, ...
    'cond_max_n', 20000, ...
    'cond_full_max_n', 1400, ...
    'cond_eigs_tol', 1e-8, ...
    'cond_eigs_maxit', 2000, ...
    'rng_seed', 20260505);
runCfg = merge_run_config_local(defaultRunCfg, runCfg);

Nc = max(1, round(runCfg.Nc));
pdeg = max(1, round(runCfg.p));
t = pdeg - 1;
refineList = unique(max(1, round(runCfg.refine_list(:).')), 'stable');
nEigenvalues = max(1, round(runCfg.n_eigenvalues));

resultRoot = fullfile(exampleDir, 'data', 'preconditioner');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

outDir = fullfile(resultRoot, 'condition_time');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cachePwRoot = fullfile(outDir, 'cache_pw');
cacheNurbsRoot = fullfile(outDir, 'cache_nurbs');
if ~exist(cachePwRoot, 'dir'), mkdir(cachePwRoot); end
if ~exist(cacheNurbsRoot, 'dir'), mkdir(cacheNurbsRoot); end

common = make_common_opts_local(runCfg, nEigenvalues, cachePwRoot, cacheNurbsRoot);

solveModes = {'none', 'purediag', 'interfaceblock'};
methodPrefixes = {'none', 'pd', 'ib'};
methodLabels = {'Unpreconditioned', 'Jacobi', 'TB-DG'};

fprintf('\n============================================================\n');
fprintf('[RUN ] Example 3 SCF preconditioner comparison\n');
fprintf('[INFO] Nc = %d, p = %d, refine = %s\n', Nc, pdeg, mat2str(refineList));
fprintf('[INFO] PW inner correction = Chebyshev, NURBS = square-fast Gauss, DG = square-fast\n');
fprintf('[INFO] output = %s\n', outDir);

rows = struct([]);
for ir = 1:numel(refineList)
    refine = refineList(ir);
    nElem = 2 ^ refine;
    hasCondCtx = false;
    row = struct();
    row.Nc = Nc;
    row.p = pdeg;
    row.refine = refine;
    row.nElem = nElem;
    row.h = 0.4 / nElem;

    for im = 1:numel(solveModes)
        opts = common;
        opts.solve_mode = solveModes{im};
        opts.rng_seed = runCfg.rng_seed + 1000 * ir + im;
        opts.outDir = fullfile(outDir, sprintf('refine_%02d', refine), solveModes{im});
        if ~exist(opts.outDir, 'dir'), mkdir(opts.outDir); end

        fprintf('\n------------------------------------------------------------\n');
        fprintf('[CASE] refine = %d, nElem = %d, method = %s\n', ...
            refine, nElem, methodLabels{im});

        runMat = fullfile(opts.outDir, 'run.mat');
        [loaded, result, meta] = load_optimized_case_local(runMat);
        if loaded
            fprintf('[CACHE] using existing optimized run.mat: %s\n', runMat);
        else
            solve_iga_pw_dg( ...
                nElem, t, Nc, nEigenvalues, opts);
            [loaded, result, meta] = load_optimized_case_local(runMat);
            assert(loaded, 'Missing optimized result data after running case: %s', runMat);
        end

        if ~hasCondCtx
            condCtx = build_condition_context_local(scriptDir, nElem, t, Nc, common);
            hasCondCtx = true;
        end
        [result, meta] = recompute_shifted_condition_local(result, meta, ...
            condCtx, opts, solveModes{im});
        update_cached_run_condition_local(runMat, result, meta);

        row = add_method_result_local(row, methodPrefixes{im}, result);
    end

    if isempty(rows)
        rows = row;
    else
        rows(end + 1, 1) = row; %#ok<AGROW>
    end
end

summary = struct2table(rows);
summaryCsv = fullfile(outDir, 'summary.csv');
writetable(summary, summaryCsv);

plotDir = fullfile(outDir, 'figures');
if ~exist(plotDir, 'dir'), mkdir(plotDir); end
delete_extra_figures_local(plotDir);
plot_preconditioner_summary_local(summary, methodPrefixes, methodLabels, plotDir);

fprintf('\n[SAVED] %s\n', summaryCsv);
fprintf('[SAVED] %s\n', fullfile(plotDir, 'time.pdf'));
fprintf('[SAVED] %s\n', fullfile(plotDir, 'condition.pdf'));
end

function opts = make_common_opts_local(runCfg, nEigenvalues, cachePwRoot, cacheNurbsRoot)
%Build common opts.
opts = struct();
opts.Example = 'Example_3';
opts.beta = runCfg.beta;
opts.n_gp = runCfg.n_gp;
opts.primme_tol = runCfg.primme_tol;
opts.primme_maxit = 1e7;
opts.primme_method = 'DEFAULT_MIN_MATVECS';
opts.primme_reportLevel = runCfg.primme_reportLevel;
opts.block_targetShift = 0.0;
opts.scf_maxit = runCfg.scf_maxit;
opts.scf_pw_grid_m = 500;
opts.scf_tol_lambda = runCfg.scf_tol_lambda;
opts.scf_mixing = 0.9;
opts.scf_track_n_eigs = max(2, nEigenvalues);
opts.scf_reuse_previous_as_v0 = true;
opts.scf_auto_damp_on_oscillation = true;
opts.eps_diag = 1e-12;
opts.iface_reg = 1e-12;
opts.rng_seed = runCfg.rng_seed;
opts.snapshot_iters = [];
opts.save_eigenvectors = false;
opts.save_nurbs = false;
opts.save_pw_index = false;
opts.save_matrices = false;
opts.save_mat = false;
opts.use_pw_cache = true;
opts.use_nurbs_cache = true;
opts.cacheRoot = cachePwRoot;
opts.cacheNurbsRoot = cacheNurbsRoot;
opts.prec_cond_max_n = runCfg.cond_max_n;
opts.inner_cheb_n = 80;
opts.pw_fft_grid_n = 128;
opts.use_square_nurbs_fast = true;
opts.use_square_dg_fast = true;
end

function [loaded, result, meta] = load_optimized_case_local(runMat)
%Load optimized case.
loaded = false;
result = struct();
meta = struct();
if ~exist(runMat, 'file')
    return;
end

S = load(runMat, 'run');
assert(isfield(S, 'run'), 'Cached file does not contain run: %s', runMat);
assert(isfield(S.run, 'result'), 'Cached run does not contain result: %s', runMat);
assert(isfield(S.run, 'meta'), 'Cached run does not contain meta: %s', runMat);

meta = S.run.meta;
result = S.run.result;
loaded = true;
end

function condCtx = build_condition_context_local(scriptDir, nElem, t, Nc, opts)
%Build condition context.
L = 4;
a = 0.2;
innerDomains = [-a, a, -a, a];

nu = 2; nv = 2;
ConPts = zeros(nu, nv, 2);
x = [-a, a];
y = [-a, a];
ConPts(:, :, 1) = [x(1) x(1); x(2) x(2)];
ConPts(:, :, 2) = [y(1) y(2); y(1) y(2)];

nurbs_original = struct();
nurbs_original.ConPts = ConPts;
nurbs_original.weights = [1 1; 1 1];
nurbs_original.pu = 1;
nurbs_original.pv = 1;
nurbs_original.knotU = [0 0 1 1];
nurbs_original.knotV = [0 0 1 1];

pu = 1 + t;
pv = 1 + t;
nurbs_refine = IGA_2D_Grid_nElem([], [], pu, pv, nElem);
nDofsNurbs = nurbs_refine.n_dofs_domains;

[kPw, nPwBasis] = build_pw_disk_local(Nc);
pwDofs = nDofsNurbs + (1:nPwBasis);
nTotal = nDofsNurbs + nPwBasis;

NVr = 1;
[kVr, nPwVr] = build_pw_disk_local(NVr);

nGp = opts.n_gp;
[~, Mnurbs] = get_nurbs_matrices_for_condition_local( ...
    nurbs_original, nurbs_refine, nElem, t, kVr, nPwVr, L, nGp, opts);
[~, Mpw] = get_pw_matrices_for_condition_local( ...
    L, Nc, innerDomains, kVr, nPwVr, opts);

M = sparse(nTotal, nTotal);
M(1:nDofsNurbs, 1:nDofsNurbs) = Mnurbs;
M(pwDofs, pwDofs) = Mpw;
M = 0.5 * (M + M');

if opts.use_square_dg_fast
    P = assemble_DG_square_interface_fast(nurbs_refine, kPw, pwDofs, L, a, nTotal);
else
    [PBottom, ~] = IGA_DG_Bottom_Edge_Assemble(nurbs_original, nurbs_refine, kPw, pwDofs, L, nTotal);
    [PTop,    ~] = IGA_DG_Top_Edge_Assemble(nurbs_original, nurbs_refine, kPw, pwDofs, L, nTotal);
    [PLeft,   ~] = IGA_DG_Left_Edge_Assemble(nurbs_original, nurbs_refine, kPw, pwDofs, L, nTotal);
    [PRight,  ~] = IGA_DG_Right_Edge_Assemble(nurbs_original, nurbs_refine, kPw, pwDofs, L, nTotal);
    P = PBottom + PTop + PLeft + PRight;
end
P = sparse(P);
gamma = find(sum(abs(P), 2) ~= 0);
eta = setdiff((1:nTotal).', gamma);

condCtx = struct();
condCtx.M = M;
condCtx.gamma = gamma;
condCtx.eta = eta;
condCtx.n_total = nTotal;
condCtx.script_dir = scriptDir;
end

function [result, meta] = recompute_shifted_condition_local(result, meta, condCtx, opts, solveMode)
%Compute shifted condition.
tau = double(opts.block_targetShift);
epsD = double(opts.eps_diag);
maxN = double(opts.prec_cond_max_n);
assert(isfield(result, 'Mat_final') && ~isempty(result.Mat_final), ...
    'Missing final operator matrix for the condition-number plot.');

A = sparse(result.Mat_final);
A = 0.5 * (A + A');
A_tau = A - tau * condCtx.M;
A_tau = 0.5 * (A_tau + A_tau');

[Pmat, ~] = build_condition_preconditioner_matrix_local(A_tau, condCtx, solveMode, epsD);
label = sprintf('%s, tau = %.6g', char(solveMode), tau);
[kappa, ~, ~] = shifted_generalized_abs_condition_local(A_tau, Pmat, maxN, label, opts);

result.cond = kappa;
end

function update_cached_run_condition_local(runMat, result, meta)
%Update cached run condition.
assert(exist(runMat, 'file') == 2, 'Missing cached run file: %s', runMat);
S = load(runMat, 'run');
assert(isfield(S, 'run'), 'Cached run file does not contain run: %s', runMat);
run = S.run;
run.result = result;
run.meta = meta;
save(runMat, 'run', '-v7.3');
end

function [Pmat, delta] = build_condition_preconditioner_matrix_local(A_tau, condCtx, solveMode, epsD)
%Build condition preconditioner matrix.
n = size(A_tau, 1);
solveMode = lower(string(solveMode));
delta = 0;

switch solveMode
    case "none"
        Pmat = speye(n);

    case "purediag"
        d = max(abs(diag(A_tau)), epsD);
        Pmat = spdiags(d, 0, n, n);

    case "interfaceblock"
        gamma = condCtx.gamma;
        eta = condCtx.eta;

        dEta = max(abs(diag(A_tau(eta, eta))), epsD);
        Ag = A_tau(gamma, gamma);
        Ag = 0.5 * (Ag + Ag');
        [delta, AgShifted] = make_block_positive_definite_local(Ag, epsD);

        Pmat = sparse(n, n);
        Pmat(eta, eta) = spdiags(dEta, 0, numel(eta), numel(eta));
        Pmat(gamma, gamma) = AgShifted;

    otherwise
        error('Unsupported solve_mode for condition number: %s', solveMode);
end

Pmat = 0.5 * (Pmat + Pmat');
end

function [delta, AgShifted] = make_block_positive_definite_local(Ag, epsD)
%Build block positive definite.
n = size(Ag, 1);
if n == 0
    delta = 0;
    AgShifted = sparse(0, 0);
    return;
end

scale = max(1, norm(Ag, 1));
floorEig = max(epsD, epsD * scale);
if n <= 3000
    lamMin = min(real(eig(full(Ag))));
else
    eigOpts = struct('issym', true, 'isreal', isreal(Ag), 'tol', 1e-8, 'maxit', 2000);
    lamMin = real(eigs(Ag, 1, 'smallestreal', eigOpts));
end

delta = max(0, floorEig - lamMin);
I = speye(n);
AgShifted = Ag + delta * I;

for attempt = 1:8
    [~, pflag] = chol(AgShifted);
    if pflag == 0
        return;
    end
    delta = max(10 * max(delta, floorEig), floorEig * 10^attempt);
    AgShifted = Ag + delta * I;
end
error('TB-DG gamma block is not numerically positive definite after delta regularization.');
end

function [kappa, sigmaMax, sigmaMin] = shifted_generalized_abs_condition_local(A, B, maxN, label, opts)
%Compute generalized abs condition.
n = size(A, 1);
assert(n <= maxN, '%s: matrix size %d exceeds max_n = %g.', label, n, maxN);

fullMaxN = double(opts.cond_full_max_n);
tol = double(opts.cond_eigs_tol);
maxit = double(opts.cond_eigs_maxit);

if n <= fullMaxN
    mu = eig(full(A), full(B), 'vector');
else
    eigOpts = struct();
    eigOpts.tol = tol;
    eigOpts.maxit = maxit;
    muMax = eigs(A, B, 1, 'largestabs', eigOpts);
    muMin = eigs(A, B, 1, 'smallestabs', eigOpts);
    mu = [muMax; muMin];
end

absMu = abs(mu(isfinite(mu)));
absMu = absMu(absMu > 0);
assert(~isempty(absMu), '%s: condition-number eigenvalues are empty.', label);
sigmaMax = max(absMu);
sigmaMin = min(absMu);
assert(sigmaMin > 0, '%s: smallest absolute eigenvalue is not positive.', label);
kappa = sigmaMax / sigmaMin;
end

function [Hnurbs, Mnurbs] = get_nurbs_matrices_for_condition_local( ...
nurbs_original, nurbs_refine, nElem, t, kVr, nPwVr, L, nGp, opts)
%Return NURBS matrices for condition.
useFast = logical(opts.use_square_nurbs_fast);
exampleName = opts.Example;
cacheFile = fullfile(opts.cacheNurbsRoot, ...
    sprintf('NURBS_%s_squarefast_%d_nElem_%03d_t_%d_ngp_%d_L_%g_NVr_%d.mat', ...
    exampleName, useFast, nElem, t, nGp, L, nPwVr));

if opts.use_nurbs_cache && exist(cacheFile, 'file')
    S = load(cacheFile, 'H_nurbs', 'M_nurbs');
    Hnurbs = S.H_nurbs;
    Mnurbs = S.M_nurbs;
    return;
end

[Hnurbs, Mnurbs, nurbsTiming] = generate_A_M_NURBS_2D( ...
    nurbs_original, nurbs_refine, kVr, nPwVr, L, nGp, exampleName, opts);
if opts.use_nurbs_cache
    if ~exist(opts.cacheNurbsRoot, 'dir'), mkdir(opts.cacheNurbsRoot); end
    H_nurbs = Hnurbs;
    M_nurbs = Mnurbs;
    save(cacheFile, 'H_nurbs', 'M_nurbs', 'nurbsTiming', '-v7.3');
end
end

function [Hpw, Mpw] = get_pw_matrices_for_condition_local( ...
L, Nc, innerDomains, kVr, nPwVr, opts)
%Return PW matrices for condition.
innerChebN = opts.inner_cheb_n;
fftGridN = opts.pw_fft_grid_n;
exampleName = opts.Example;
cacheFile = fullfile(opts.cacheRoot, ...
    sprintf('PW_chebfast_%s_fft_%d_cheb_%d_L_%g_Nc_%d_NVr_%d.mat', ...
    exampleName, fftGridN, innerChebN, L, Nc, nPwVr));

if opts.use_pw_cache && exist(cacheFile, 'file')
    S = load(cacheFile, 'H_pw', 'M_pw');
    Hpw = S.H_pw;
    Mpw = S.M_pw;
    return;
end

[Hpw, Mpw, ~, pwTiming] = generate_A_M_PW_2D(L, Nc, innerDomains, kVr, nPwVr, opts);
if opts.use_pw_cache
    if ~exist(opts.cacheRoot, 'dir'), mkdir(opts.cacheRoot); end
    H_pw = Hpw;
    M_pw = Mpw;
    save(cacheFile, 'H_pw', 'M_pw', 'pwTiming', '-v7.3');
end
end

function [kList, nBasis] = build_pw_disk_local(Nc)
%Build PW disk.
N = floor(Nc);
kList = zeros((2 * N + 1)^2, 2);
nBasis = 0;
for i = -N:N
    m = floor(sqrt(N^2 - i^2));
    for j = -m:m
        nBasis = nBasis + 1;
        kList(nBasis, :) = [i, j];
    end
end
kList = kList(1:nBasis, :);
end

function row = add_method_result_local(row, prefix, result)
%Store only values used by the manuscript plots.
row.([prefix '_solver_time_s']) = result.time_total;
row.([prefix '_cond']) = result.cond;
end

function plot_preconditioner_summary_local(summary, prefixes, labels, plotDir)
%Plot preconditioner summary.
cfg = default_style_local();
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

[hvals, order] = sort(summary.h, 'ascend');
solverTotal = extract_metric_matrix_local(summary, prefixes, '_solver_time_s', order);
condv = extract_metric_matrix_local(summary, prefixes, '_cond', order);

fig1 = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);
ax = axes(fig1, 'Position', [cfg.layout.left, cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top]);
plot_metric_panel_local(ax, hvals, solverTotal, labels, cfg, ...
    '$h$', '$\mathrm{Solver\ total\ time\ (s)}$');
legend(ax, labels, ...
    'Location', cfg.legend.location, ...
    'Box', cfg.legend.box, ...
    'FontSize', cfg.legend.fontSize, ...
    'Interpreter', 'latex');
export_figure_local(fig1, fullfile(plotDir, 'time'), cfg);
close(fig1);

fig2 = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);
ax = axes(fig2, 'Position', [cfg.layout.left, cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top]);
plot_metric_panel_local(ax, hvals, condv, labels, cfg, ...
    '$h$', '$\mathrm{Condition\ Number}$');
ylim(ax, [1, 1e10]);
set(ax, ...
    'YTick', 10 .^ [0 2 4 6 8 10], ...
    'YTickLabel', {'$10^0$', '$10^2$', '$10^4$', '$10^6$', '$10^8$', '$10^{10}$'});
legend(ax, labels, ...
    'Location', cfg.legend.location, ...
    'Box', cfg.legend.box, ...
    'FontSize', cfg.legend.fontSize, ...
    'Interpreter', 'latex');
export_figure_local(fig2, fullfile(plotDir, 'condition'), cfg);
close(fig2);
end

function M = extract_metric_matrix_local(summary, prefixes, suffix, order)
%Extract metric matrix.
M = zeros(numel(order), numel(prefixes));
for k = 1:numel(prefixes)
    M(:, k) = summary.([prefixes{k} suffix])(order);
end
assert(all(isfinite(M(:)) & M(:) > 0), 'Summary data must be finite and positive.');
end

function plot_metric_panel_local(ax, x, Y, labels, cfg, xlabelStr, ylabelStr)
%Plot metric panel.
X = repmat(x(:), 1, numel(labels));
plot_xy_metric_panel_local(ax, X, Y, labels, cfg, xlabelStr, ylabelStr);
end

function plot_xy_metric_panel_local(ax, X, Y, labels, cfg, xlabelStr, ylabelStr)
%Plot xy metric panel.
hold(ax, 'on');
for k = 1:numel(labels)
    xk = X(:, k);
    yk = Y(:, k);
    assert(all(isfinite(xk) & isfinite(yk) & xk > 0 & yk > 0), ...
        'Plot data must be finite and positive.');
    plot(ax, xk, yk, '-', ...
        'LineWidth', cfg.line.width, ...
        'Color', cfg.line.colors(k, :), ...
        'Marker', cfg.line.markers{k}, ...
        'MarkerSize', cfg.line.markerSize, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', cfg.line.colors(k, :));
end
set(ax, 'XScale', 'log', 'YScale', 'log');
apply_axis_padding_local(ax, X, Y, cfg);
set_axes_style_local(ax, cfg);
xlabel(ax, xlabelStr, 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
ylabel(ax, ylabelStr, 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
end

function apply_axis_padding_local(ax, x, Y, cfg)
%Add padding to axis limits.
x = x(:);
y = Y(:);
assert(all(isfinite(x) & x > 0), 'X data for axis padding must be finite and positive.');
assert(all(isfinite(y) & y > 0), 'Y data for axis padding must be finite and positive.');
lx = log10([min(x), max(x)]);
dx = max(diff(lx), eps);
xlim(ax, 10 .^ [lx(1) - cfg.axes.xpad * dx, lx(2) + cfg.axes.xpad * dx]);
ly = log10([min(y), max(y)]);
dy = max(diff(ly), eps);
ylim(ax, 10 .^ [ly(1) - cfg.axes.ypad * dy, ly(2) + cfg.axes.ypad * dy]);
end

function set_axes_style_local(ax, cfg)
%Apply axes style settings.
set(ax, ...
    'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'Box', 'on', ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick, ...
    'XMinorGrid', 'off', ...
    'YMinorGrid', 'off');
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
if isprop(ax, 'XRuler')
    ax.XRuler.MinorTick = 'off';
end
if isprop(ax, 'YRuler')
    ax.YRuler.MinorTick = 'off';
end
grid(ax, 'off');
end

function export_figure_local(fig, baseName, cfg)
%Export figure.
set(fig, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0 0 cfg.fig.width cfg.fig.height], ...
    'PaperSize', [cfg.fig.width cfg.fig.height], ...
    'InvertHardcopy', 'off');
exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
end

function cfg = default_style_local()
%Return plotting style values.
cfg = struct();
cfg.fig.width = 4.8;
cfg.fig.height = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor = 'w';
cfg.layout.left = 0.14;
cfg.layout.right = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top = 0.08;
cfg.axes.fontSize = 10;
cfg.axes.lineWidth = 1.0;
cfg.axes.tickDir = 'out';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';
cfg.axes.labelSize = 12;
cfg.axes.xpad = 0.06;
cfg.axes.ypad = 0.08;
cfg.legend.location = 'northeast';
cfg.legend.box = 'off';
cfg.legend.fontSize = 11;
cfg.line.colors = [ ...
    223 122 094;
    060 064 091;
    130 178 154] / 255;
cfg.line.markers = {'o', 's', '^'};
cfg.line.width = 1.8;
cfg.line.markerSize = 8;
end

function delete_extra_figures_local(plotDir)
%Close figures that are not needed.
keep = {'time.pdf', 'condition.pdf'};
files = dir(fullfile(plotDir, '*'));
for k = 1:numel(files)
    if files(k).isdir
        continue;
    end
    if ~ismember(files(k).name, keep)
        delete(fullfile(files(k).folder, files(k).name));
    end
end
end

function cfg = merge_run_config_local(defaultCfg, userCfg)
%Merge run config.
cfg = defaultCfg;
assert(isstruct(userCfg), 'Run configuration must be a structure.');
names = fieldnames(userCfg);
for k = 1:numel(names)
    cfg.(names{k}) = userCfg.(names{k});
end
end
