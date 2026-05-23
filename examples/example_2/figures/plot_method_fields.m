function plot_method_fields()
%Plot density and potential fields.

clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);

activate_example_workflow('method_comparison', ...
    {'nurbs', 'iga', 'operators', 'core'});

cfgPW = struct('Nc', 20);
cfgIGA = struct('pdeg', 1, 'refine', 6);
cfgHybrid = struct('Nc', 10, 'pdeg', 1, 'refine', 4, 't', 0);
cfgRef = struct('Nc', 45, 'pdeg', 3, 'refine', 8, 't', 2);

outDir = fullfile(repoRoot, 'data', 'method_fields');
refDir = fullfile(outDir, 'reference');
pwDir = fullfile(outDir, 'pw');
igaDir = fullfile(outDir, 'iga');
igapwDir = fullfile(outDir, 'iga_pw');
fieldsDir = fullfile(outDir, 'fields');
ensure_dir_local(outDir);
ensure_dir_local(refDir);
ensure_dir_local(pwDir);
ensure_dir_local(igaDir);
ensure_dir_local(igapwDir);
ensure_dir_local(fieldsDir);

common = make_common_opts_local(scriptDir);

fprintf('\n============================================================\n');
fprintf('[RUN ] Example 2 PW / IGA / IGA-PW u1 comparison\n');

refRun = load_or_run_reference_local(refDir, cfgRef, common);
pwRun = load_or_run_pw_local(pwDir, cfgPW, common);
igaRun = load_or_run_iga_local(igaDir, cfgIGA, common);
igapwRun = load_or_run_hybrid_local(igapwDir, cfgHybrid, common);

gridN = 401;
[Xg, Yg] = meshgrid(linspace(-2, 2, gridN), linspace(-2, 2, gridN));
dxg = Xg(1, 2) - Xg(1, 1);

Uref_raw = evaluate_run_on_grid_local(refRun, 'igapw', Xg, Yg);
[Uref, refNormMeta] = normalize_and_align_local(Uref_raw, [], dxg);

Upw_raw = evaluate_run_on_grid_local(pwRun, 'pw', Xg, Yg);
[Upw, pwNormMeta] = normalize_and_align_local(Upw_raw, Uref, dxg);

Uiga_raw = evaluate_run_on_grid_local(igaRun, 'iga', Xg, Yg);
[Uiga, igaNormMeta] = normalize_and_align_local(Uiga_raw, Uref, dxg);

Uigapw_raw = evaluate_run_on_grid_local(igapwRun, 'igapw', Xg, Yg);
[Uigapw, igapwNormMeta] = normalize_and_align_local(Uigapw_raw, Uref, dxg);

U = {Upw; Uiga; Uigapw; Uref};
Err = {abs(Upw - Uref); abs(Uiga - Uref); abs(Uigapw - Uref)};
methodNames = {'PW', 'IGA', 'IGA-PW'};
plotNames = {'PW', 'IGA', 'IGA-PW', 'Reference'};

summaryTable = build_summary_table_local(pwRun, igaRun, igapwRun, refRun, cfgPW, cfgIGA, cfgHybrid, cfgRef, Err, dxg);
writetable(summaryTable, fullfile(outDir, 'summary.csv'));

save(fullfile(fieldsDir, 'fields.mat'), ...
    'Xg', 'Yg', 'U', 'Err', 'refNormMeta', 'pwNormMeta', 'igaNormMeta', 'igapwNormMeta', ...
    'pwRun', 'igaRun', 'igapwRun', 'refRun', '-v7.3');

cleanup_legacy_outputs_local(outDir);

example1RefCmap = example1_reference_colormap_local(256);
fieldClimU1 = shared_field_clim_local(U(1:3));
errorClimU1 = shared_error_clim_local(Err);

export_split_field_set_local(Xg, Yg, U, plotNames, outDir, 'u1', example1RefCmap, fieldClimU1);
export_split_error_set_local(Xg, Yg, Err, methodNames, outDir, 'u1', example1RefCmap, errorClimU1);

write_caption_file_local(outDir, cfgPW, cfgIGA, cfgHybrid, cfgRef, summaryTable);
write_report_local(outDir, summaryTable, cfgPW, cfgIGA, cfgHybrid, cfgRef);

fprintf('[DONE] Example 2 u1 comparison outputs saved to:\n%s\n', outDir);
fprintf('============================================================\n\n');

end

function common = make_common_opts_local(scriptDir)
%Build common opts.
exampleDir = fileparts(scriptDir);
common = struct();
common.Example = 'Example_2';
common.beta = 20;
common.n_gp = 10;
common.inner_cheb_n = 48;
common.pw_fft_grid_n = 256;
common.primme_tol = 1e-12;
common.primme_maxit = 2e8;
common.primme_method = 'DEFAULT_MIN_MATVECS';
common.primme_reportLevel = 0;
common.n_eigenvalues = 1;
common.save_eigenvectors = true;
common.save_nurbs = true;
common.save_pw_index = true;
common.save_matrices = false;
common.save_mat = false;
common.use_pw_cache = true;
common.cacheRoot = fullfile(exampleDir, 'data', 'method_fields', 'cache_pw_shared');
common.cacheNurbsRoot = fullfile(exampleDir, 'data', 'method_fields', 'cache_nurbs_shared');
common.cacheInterfaceRoot = fullfile(exampleDir, 'data', 'method_fields', 'cache_interface_shared');
common.potential_Nc_pw = 2;
common.potential_Nc_iga = 2;
common.pw_potential_grid_m = 1024;
common.eps_diag = 1e-12;
common.iface_reg = 1e-12;
ensure_dir_local(common.cacheRoot);
ensure_dir_local(common.cacheNurbsRoot);
ensure_dir_local(common.cacheInterfaceRoot);
end

function run = load_or_run_reference_local(outDir, cfg, common)
%Load or run reference.
fprintf('[CASE] reference IGA-PW (Nc=%d, p=%d, refine=%d, h-reference)\n', cfg.Nc, cfg.pdeg, cfg.refine);
runFile = fullfile(outDir, 'run.mat');

exampleDir = fileparts(fileparts(mfilename('fullpath')));
latestRunFile = fullfile(exampleDir, 'data', 'result', 'Example_2', ...
    sprintf('Nc_%d', cfg.Nc), sprintf('p_%d', cfg.pdeg), ...
    sprintf('refine_%02d', cfg.refine), 'run.mat');
assert(exist(latestRunFile, 'file') == 2, 'Missing h-reference run: %s', latestRunFile);
S = load(latestRunFile, 'run');
assert(isfield(S, 'run'), 'Missing run structure in h-reference: %s', latestRunFile);
run = S.run;
assert(isstruct(run) && isfield(run, 'meta') && isstruct(run.meta), ...
    'Reference run metadata is missing: %s', latestRunFile);
assert(isfield(run.meta, 'Nc') && isequal(double(run.meta.Nc), double(cfg.Nc)), ...
    'Reference Nc does not match requested Nc=%d.', cfg.Nc);
assert(isfield(run.meta, 'Refinement') && isequal(double(run.meta.Refinement), double(cfg.refine)), ...
    'Reference refine does not match requested refine=%d.', cfg.refine);
assert(isfield(run.meta, 'pu') && isequal(double(run.meta.pu), double(cfg.pdeg)), ...
    'Reference p does not match requested p=%d.', cfg.pdeg);
assert(isfield(run, 'n_dofs_1') && isfield(run, 'n_dofs_2') ...
    && isfield(run, 'n_dofs_nurbs') && isfield(run, 'n_pw_basis') ...
    && isfield(run, 'k_pw') && isfield(run, 'nurbs_refine_1') ...
    && isfield(run, 'nurbs_refine_2') && isfield(run, 'uh') ...
    && ~isempty(run.uh) && isfield(run, 'lambda') && numel(run.lambda) >= 1, ...
    'Reference run lacks the field data needed for method comparison.');
save(runFile, 'run', '-v7.3');
end

function run = load_or_run_pw_local(outDir, cfg, common)
%Load or run PW.
fprintf('[CASE] PW-only (Nc=%d)\n', cfg.Nc);
runFile = fullfile(outDir, 'run.mat');
if exist(runFile, 'file')
    S = load(runFile, 'run');
    if is_single_run_usable_local(S.run, 'pw', 'purediag', cfg.Nc, 0, false)
        run = S.run;
        return;
    end
end
opts = common;
opts.outDir = outDir;
opts.solve_mode = 'purediag';
solve_full_domain('pw', struct('Nc', cfg.Nc), opts);
S = load(runFile, 'run');
run = S.run;
end

function run = load_or_run_iga_local(outDir, cfg, common)
%Load or run IGA.
fprintf('[CASE] IGA-only (p=%d, refine=%d)\n', cfg.pdeg, cfg.refine);
runFile = fullfile(outDir, 'run.mat');
if exist(runFile, 'file')
    S = load(runFile, 'run');
    if is_single_run_usable_local(S.run, 'iga', 'purediag', cfg.refine, cfg.pdeg, true)
        run = S.run;
        return;
    end
end
opts = common;
opts.outDir = outDir;
opts.solve_mode = 'purediag';
opts.periodic_bc_iga = true;
solve_full_domain('iga', struct('pdeg', cfg.pdeg, 'refine', cfg.refine), opts);
S = load(runFile, 'run');
run = S.run;
end

function run = load_or_run_hybrid_local(outDir, cfg, common)
%Load or run hybrid.
fprintf('[CASE] IGA-PW (Nc=%d, p=%d, refine=%d, InterfaceBlock)\n', cfg.Nc, cfg.pdeg, cfg.refine);
runFile = fullfile(outDir, 'run.mat');
if exist(runFile, 'file')
    S = load(runFile, 'run');
    if is_hybrid_run_usable_local(S.run, 'interfaceblock', cfg.Nc, cfg.refine, cfg.pdeg)
        run = S.run;
        return;
    end
end
opts = common;
opts.outDir = outDir;
opts.solve_mode = 'interfaceblock';
solve_iga_pw_dg(cfg.refine, cfg.t, cfg.Nc, 1, opts);
S = load(runFile, 'run');
run = S.run;
end

function summaryTable = build_summary_table_local(pwRun, igaRun, igapwRun, refRun, cfgPW, cfgIGA, cfgHybrid, cfgRef, Err, dxg)
%Build summary table.
method = {'PW'; 'IGA'; 'IGA-PW'; 'Reference'};
Nc = [cfgPW.Nc; 0; cfgHybrid.Nc; cfgRef.Nc];
pdeg = [0; cfgIGA.pdeg; cfgHybrid.pdeg; cfgRef.pdeg];
refine = [0; cfgIGA.refine; cfgHybrid.refine; cfgRef.refine];
h = [0; igaRun.meta.h; igapwRun.meta.h; refRun.meta.h];
lambda1 = [pwRun.lambda(1); igaRun.lambda(1); igapwRun.lambda(1); refRun.lambda(1)];
abs_err_lambda1 = abs(lambda1 - refRun.lambda(1));
l2_error_u1 = [l2_error_from_grid_local(Err{1}, dxg); l2_error_from_grid_local(Err{2}, dxg); l2_error_from_grid_local(Err{3}, dxg); 0];
summaryTable = table(method, Nc, pdeg, refine, h, lambda1, abs_err_lambda1, l2_error_u1);
end

function val = l2_error_from_grid_local(fieldVals, dx)
%Compute error from grid.
val = sqrt(sum(abs(fieldVals(:)).^2) * dx * dx);
end

function [Uout, meta] = normalize_and_align_local(Uin, Uref, dx)
%Normalize and align.
u = Uin;
norm_before = sqrt(max(sum(abs(u(:)).^2) * dx * dx, eps));
u = u ./ norm_before;

phase_applied = 0;
if isempty(Uref)
    [~, idx] = max(abs(u(:)));
    phase_applied = angle(u(idx));
    u = u * exp(-1i * phase_applied);
    if real(u(idx)) < 0
        u = -u;
        phase_applied = phase_applied + pi;
    end
else
    innerVal = sum(u(:) .* conj(Uref(:))) * dx * dx;
    if abs(innerVal) > 0
        phase_applied = angle(innerVal);
        u = u * exp(-1i * phase_applied);
    end
    innerVal2 = real(sum(u(:) .* conj(Uref(:))) * dx * dx);
    if innerVal2 < 0
        u = -u;
        phase_applied = phase_applied + pi;
    end
end

Uout = u;
meta = struct('norm_before', norm_before, 'phase_applied', phase_applied);
end

function Ugrid = evaluate_run_on_grid_local(run, kind, Xg, Yg)
%Evaluate a saved solution on a plotting grid.
coeff = run.uh(:, 1);
L = 4;

switch lower(kind)
    case 'pw'
        Ugrid = reshape(pw_eval_val_chunked_local(coeff, run.k_pw, Xg(:), Yg(:), L, 20000), size(Xg));

    case 'iga'
        Ugrid = reshape(iga_eval_on_one_patch_local(run.nurbs_refine, coeff, Xg(:), Yg(:), 0.0, 0.0, 2.0), size(Xg));

    case 'igapw'
        patchCenters = [-1, 0; 1, 0];
        a = 0.2;
        n1 = run.n_dofs_1;
        n2 = run.n_dofs_2;
        coeff1 = coeff(1:n1);
        coeff2 = coeff(n1 + (1:n2));
        coeffPw = coeff(run.n_dofs_nurbs + (1:run.n_pw_basis));

        x = Xg(:);
        y = Yg(:);
        mask1 = (x >= patchCenters(1,1)-a) & (x <= patchCenters(1,1)+a) ...
            & (y >= patchCenters(1,2)-a) & (y <= patchCenters(1,2)+a);
        mask2 = (x >= patchCenters(2,1)-a) & (x <= patchCenters(2,1)+a) ...
            & (y >= patchCenters(2,2)-a) & (y <= patchCenters(2,2)+a);
        maskOuter = ~(mask1 | mask2);

        val = zeros(size(x));
        if any(maskOuter)
            val(maskOuter) = pw_eval_val_chunked_local(coeffPw, run.k_pw, x(maskOuter), y(maskOuter), L, 20000);
        end
        if any(mask1)
            val(mask1) = iga_eval_on_one_patch_local(run.nurbs_refine_1, coeff1, x(mask1), y(mask1), patchCenters(1,1), patchCenters(1,2), a);
        end
        if any(mask2)
            val(mask2) = iga_eval_on_one_patch_local(run.nurbs_refine_2, coeff2, x(mask2), y(mask2), patchCenters(2,1), patchCenters(2,2), a);
        end
        Ugrid = reshape(val, size(Xg));

    otherwise
        error('Unknown kind: %s', kind);
end
end

function val = pw_eval_val_chunked_local(coeff, p_vec, X, Y, L, chunkSize)
%Evaluate field values in chunks.
nPts = numel(X);
val = zeros(nPts, 1);
for i1 = 1:chunkSize:nPts
    i2 = min(i1 + chunkSize - 1, nPts);
    F = [X(i1:i2).'; Y(i1:i2).'];
    expo = exp((1i * 2 * pi / L) * (p_vec * F));
    val(i1:i2) = ((coeff.' * expo) / L).';
end
end

function val = iga_eval_on_one_patch_local(nurbs, coeff, X, Y, xc, yc, a)
%Compute eval on one patch.
pu = nurbs.pu;
pv = nurbs.pv;
U = nurbs.Ubar(:).';
V = nurbs.Vbar(:).';

mU = length(U) - pu - 1;
nV = length(V) - pv - 1;

val = zeros(numel(X), 1);
for k = 1:numel(X)
    u = (X(k) - (xc - a)) / (2 * a);
    v = (Y(k) - (yc - a)) / (2 * a);
    u = max(0, min(1, u));
    v = max(0, min(1, v));

    spanU = findspan_local(mU - 1, pu, u, U);
    spanV = findspan_local(nV - 1, pv, v, V);

    Nu = bspline_basis_local(U, pu, u, spanU);
    Nv = bspline_basis_local(V, pv, v, spanV);

    s = 0.0 + 0.0i;
    for j1 = (spanV - pv):spanV
        lv = j1 - (spanV - pv) + 1;
        for i1 = (spanU - pu):spanU
            lu = i1 - (spanU - pu) + 1;
            row = i1 + (j1 - 1) * mU;
            s = s + coeff(row) * Nu(lu) * Nv(lv);
        end
    end
    val(k) = s;
end
end

function N = bspline_basis_local(U, p, u, span)
%Evaluate local B-spline basis values.
ndu = zeros(p + 1, p + 1);
left = zeros(1, p + 1);
right = zeros(1, p + 1);

ndu(1,1) = 1.0;
for j = 1:p
    left(j + 1) = u - U(span + 1 - j);
    right(j + 1) = U(span + j) - u;
    saved = 0.0;
    for r = 0:(j - 1)
        ndu(j + 1, r + 1) = right(r + 2) + left(j - r + 1);
        temp = ndu(r + 1, j) / ndu(j + 1, r + 1);
        ndu(r + 1, j + 1) = saved + right(r + 2) * temp;
        saved = left(j - r + 1) * temp;
    end
    ndu(j + 1, j + 1) = saved;
end
N = ndu(1:p + 1, p + 1).';
end

function span = findspan_local(n, p, u, U)
%Locate an index or object used by the computation.
if u >= U(n + 2)
    span = n + 1;
    return;
end
if u <= U(p + 1)
    span = p + 1;
    return;
end

low = p + 1;
high = n + 2;
mid = floor((low + high) / 2);
while (u < U(mid) || u >= U(mid + 1))
    if u < U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high) / 2);
end
span = mid;
end

function export_split_field_set_local(Xg, Yg, fields, labels, outDir, eigTag, cmap, climVals)
%Export split field set.
for i = 1:numel(labels)
    label = labels{i};
    if strcmpi(label, 'Reference')
        continue;
    end
    keepColorbar = strcmpi(label, 'IGA-PW');
    baseName = fullfile(outDir, build_split_field_name_local(label, eigTag));
    export_single_field_panel_local(Xg, Yg, real(fields{i}), baseName, cmap, climVals, keepColorbar);
end
end

function export_split_error_set_local(Xg, Yg, fields, labels, outDir, eigTag, cmap, climVals)
%Export split error set.
for i = 1:numel(labels)
    label = labels{i};
    keepColorbar = strcmpi(label, 'IGA-PW');
    baseName = fullfile(outDir, build_split_error_name_local(label, eigTag));
    export_single_error_panel_local(Xg, Yg, fields{i}, baseName, cmap, climVals, keepColorbar);
end
end

function export_single_field_panel_local(Xg, Yg, fieldVals, baseName, cmap, climVals, keepColorbar)
%Export single field panel.
if keepColorbar
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 3.30, 2.55]);
    axPos = [0.01, 0.01, 0.758, 0.98];
    cbPos = [0.798, 0.11, 0.05, 0.78];
else
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 2.55, 2.55]);
    axPos = [0.01, 0.01, 0.98, 0.98];
    cbPos = [];
end

ax = axes(fig, 'Position', axPos);
imagesc(ax, Xg(1, :), Yg(:, 1), fieldVals);
axis(ax, 'image');
set(ax, 'YDir', 'normal');
colormap(ax, cmap);
caxis(ax, climVals);
style_clean_image_axes_local(ax);

if keepColorbar
    draw_clean_linear_colorbar_local(fig, cbPos, cmap, climVals);
end

export_figure_raster_local(fig, baseName, 600);
close(fig);
end

function export_single_error_panel_local(Xg, Yg, fieldVals, baseName, cmap, climVals, keepColorbar)
%Export single error panel.
logMin = log10(max(climVals(1), eps));
logMax = log10(climVals(2));
errFloor = 10^logMin;
levels = linspace(logMin, logMax, 20);

if keepColorbar
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 3.30, 2.55]);
    axPos = [0.01, 0.01, 0.758, 0.98];
    cbPos = [0.798, 0.11, 0.05, 0.78];
else
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 2.55, 2.55]);
    axPos = [0.01, 0.01, 0.98, 0.98];
    cbPos = [];
end

ax = axes(fig, 'Position', axPos);
Flog = log10(max(fieldVals, errFloor));
contourf(ax, Xg, Yg, Flog, levels, 'LineStyle', 'none');
axis(ax, 'image');
set(ax, 'YDir', 'normal');
colormap(ax, cmap);
caxis(ax, [logMin, logMax]);
style_clean_image_axes_local(ax);

if keepColorbar
    draw_clean_log_colorbar_local(fig, cbPos, cmap, logMin, logMax);
end

export_figure_raster_local(fig, baseName, 600);
close(fig);
end

function climVals = shared_field_clim_local(fields)
%Compute field clim.
vals = [];
for i = 1:numel(fields)
    vals = [vals; real(fields{i}(:))]; %#ok<AGROW>
end
vals = vals(isfinite(vals));
if isempty(vals)
    climVals = [-1, 1];
    return;
end
cmin = min(vals);
cmax = max(vals);
if cmax <= cmin
    pad = max(abs(cmax), 1) * 1e-3;
    climVals = [cmin - pad, cmax + pad];
elseif cmin < 0 && cmax > 0
    amax = max(abs([cmin, cmax]));
    climVals = [-amax, amax];
else
    climVals = [cmin, cmax];
end
end

function climVals = shared_error_clim_local(fields)
%Compute error clim.
vals = [];
for i = 1:numel(fields)
    vals = [vals; abs(fields{i}(:))]; %#ok<AGROW>
end
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals)
    climVals = [eps, 1];
    return;
end
cmax = max(vals);
climVals = [max(cmax * 1e-5, eps), cmax];
end

function name = build_split_field_name_local(label, eigTag)
%Build split field name.
switch string(label)
    case "IGA"
        prefix = 'IGA_eigenfunction';
    case "IGA-PW"
        prefix = 'IGA-PW_eigenfunction';
    case "PW"
        prefix = 'PW_eigenfunction';
    otherwise
        prefix = sprintf('%s_eigenfunction', char(label));
end
name = sprintf('%s_%s', prefix, eigTag);
end

function name = build_split_error_name_local(label, eigTag)
%Build split error name.
switch string(label)
    case "IGA"
        prefix = 'IGA_error';
    case "IGA-PW"
        prefix = 'IGA-PW_error';
    case "PW"
        prefix = 'PW_error';
    otherwise
        prefix = sprintf('%s_error', char(label));
end
name = sprintf('%s_%s', prefix, eigTag);
end

function style_clean_image_axes_local(ax)
%Apply clean image axes.
set(ax, 'XTick', [], 'YTick', [], 'Box', 'off', 'LineWidth', 0.5);
axis(ax, 'off');
end

function draw_clean_linear_colorbar_local(fig, pos, cmap, climVals)
%Draw clean linear colorbar.
cbax = axes(fig, 'Position', pos);
grad = repmat(linspace(climVals(1), climVals(2), 512).', 1, 2);
imagesc(cbax, [0, 1], [climVals(1), climVals(2)], grad);
set(cbax, 'YDir', 'normal');
colormap(cbax, cmap);
set(cbax, 'XTick', [], 'YAxisLocation', 'right', 'Box', 'off', ...
    'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 0.5);
ylim(cbax, [climVals(1), climVals(2)]);
xlim(cbax, [0, 1]);
yticks(cbax, default_colorbar_ticks_local(climVals(1), climVals(2)));
end

function draw_clean_log_colorbar_local(fig, pos, cmap, logMin, logMax)
%Draw clean log colorbar.
cbax = axes(fig, 'Position', pos);
grad = repmat(linspace(logMin, logMax, 512).', 1, 2);
imagesc(cbax, [0, 1], [logMin, logMax], grad);
set(cbax, 'YDir', 'normal');
colormap(cbax, cmap);
set(cbax, 'XTick', [], 'YAxisLocation', 'right', 'Box', 'off', ...
    'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 0.5);
ylim(cbax, [logMin, logMax]);
xlim(cbax, [0, 1]);
[ticks, tickLabels] = log_colorbar_ticks_local(logMin, logMax);
yticks(cbax, ticks);
yticklabels(cbax, tickLabels);
end

function tickVals = default_colorbar_ticks_local(cmin, cmax)
%Compute colorbar ticks.
tmpFig = figure('Visible', 'off', 'Color', 'w');
tmpAx = axes(tmpFig);
imagesc(tmpAx, [0 1], [cmin cmax], [cmin cmin; cmax cmax]);
set(tmpAx, 'YDir', 'normal');
colormap(tmpAx, example1_reference_colormap_local(256));
caxis(tmpAx, [cmin, cmax]);
tmpCb = colorbar(tmpAx);
drawnow;
tickVals = tmpCb.Ticks;
if isempty(tickVals)
    tickVals = linspace(cmin, cmax, 5);
end
close(tmpFig);
end

function [ticks, tickLabels] = log_colorbar_ticks_local(logMin, logMax)
%Compute colorbar ticks.
eMin = ceil(logMin - 1e-10);
eMax = floor(logMax + 1e-10);
ticks = eMin:eMax;
if isempty(ticks)
    ticks = linspace(logMin, logMax, 3);
    tickLabels = arrayfun(@(e) sprintf('10^{%.1f}', e), ticks, 'UniformOutput', false);
    return;
end
if numel(ticks) > 6
    step = ceil(numel(ticks) / 6);
    ticks = ticks(1:step:end);
end
tickLabels = arrayfun(@(e) sprintf('10^{%d}', e), ticks, 'UniformOutput', false);
end

function export_figure_raster_local(fig, baseName, dpi)
%Export figure raster.
axList = findall(fig, 'Type', 'axes');
for iax = 1:numel(axList)
    ax = axList(iax);
    if isprop(ax, 'Toolbar')
        ax.Toolbar.Visible = 'off';
    end
end
drawnow;
set(fig, 'PaperPositionMode', 'auto');
exportgraphics(fig, [baseName '.png'], 'Resolution', dpi, 'BackgroundColor', 'white');
exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'image', 'Resolution', dpi, 'BackgroundColor', 'white');
end

function write_caption_file_local(outDir, cfgPW, cfgIGA, cfgHybrid, cfgRef, summaryTable)
%Write caption file.
fid = fopen(fullfile(outDir, 'captions_example2_pw_iga_igapw.tex'), 'w');
assert(fid >= 0, 'Failed to open caption file.');

l2u1 = summaryTable.l2_error_u1(1:3);

fprintf(fid, '%% Figure 1\n');
fprintf(fid, ['\\caption{(Example 2) Split exports of the first eigenfunction for PW ($Nc = %d$), IGA ($p = %d$, $\\mathrm{refine} = %d$), ' ...
    'and IGA-PW ($Nc = %d$, $p = %d$, $\\mathrm{refine} = %d$). The fields are evaluated on the common Cartesian grid over $[-2,2]^2$, ' ...
    '$L^2$-normalized, phase-aligned with the IGA-PW reference solution ($Nc = %d$, $p = %d$, $\\mathrm{refine} = %d$), ' ...
    'and plotted with the Example 1 reference eigenfunction colormap using a common color scale for PW, IGA, and IGA-PW.}\n\n'], ...
    cfgPW.Nc, cfgIGA.pdeg, cfgIGA.refine, cfgHybrid.Nc, cfgHybrid.pdeg, cfgHybrid.refine, cfgRef.Nc, cfgRef.pdeg, cfgRef.refine);

fprintf(fid, '%% Figure 2\n');
fprintf(fid, ['\\caption{Absolute error of the first eigenfunction: plane wave, IGA, and IGA-PW. ' ...
    'The corresponding $L^2$ errors are $%.4e$, $%.4e$, and $%.4e$, respectively, ' ...
    'computed against the IGA-PW reference solution with $Nc = %d$, $p = %d$, and $\\mathrm{refine} = %d$. ' ...
    'Each split export shows $\\log_{10}|u_1-u_{1,\\mathrm{ref}}|$ using a common logarithmic color scale.}\n'], ...
    l2u1(1), l2u1(2), l2u1(3), cfgRef.Nc, cfgRef.pdeg, cfgRef.refine);

fclose(fid);
end

function write_report_local(outDir, summaryTable, cfgPW, cfgIGA, cfgHybrid, cfgRef)
%Write a text report.
fid = fopen(fullfile(outDir, 'report_example2_pw_iga_igapw.txt'), 'w');
assert(fid >= 0, 'Failed to open report file.');

fprintf(fid, 'Example 2 PW / IGA / IGA-PW u1 comparison\n');
fprintf(fid, '==============================================\n\n');
fprintf(fid, 'Comparison cases\n');
fprintf(fid, '- PW: Nc = %d, periodic boundary condition, solve_mode = PureDiag-Jacobi\n', cfgPW.Nc);
fprintf(fid, '- IGA: p = %d, refine = %d, periodic boundary condition, solve_mode = PureDiag-Jacobi\n', cfgIGA.pdeg, cfgIGA.refine);
fprintf(fid, '- IGA-PW: Nc = %d, p = %d, refine = %d, solve_mode = InterfaceBlock\n', cfgHybrid.Nc, cfgHybrid.pdeg, cfgHybrid.refine);
fprintf(fid, '- Reference IGA-PW: Nc = %d, p = %d, refine = %d, solve_mode = InterfaceBlock\n\n', cfgRef.Nc, cfgRef.pdeg, cfgRef.refine);
fprintf(fid, 'Saved figures\n');
fprintf(fid, '- IGA_eigenfunction_u1.pdf/png\n');
fprintf(fid, '- IGA-PW_eigenfunction_u1.pdf/png\n');
fprintf(fid, '- PW_eigenfunction_u1.pdf/png\n');
fprintf(fid, '- IGA_error_u1.pdf/png\n');
fprintf(fid, '- IGA-PW_error_u1.pdf/png\n');
fprintf(fid, '- PW_error_u1.pdf/png\n\n');

fprintf(fid, 'Eigenvalue and u1 error summary\n');
fprintf(fid, '%-10s %-8s %-8s %-8s %-16s %-16s %-16s %-16s\n', ...
    'Method', 'Nc', 'p', 'refine', 'h', 'lambda1', '|lambda1-lref|', 'L2err-u1');
fprintf(fid, '%s\n', repmat('-', 1, 112));
for i = 1:height(summaryTable)
    fprintf(fid, '%-10s %-8s %-8s %-8s %-16s %-16.10f %-16.6e %-16.6e\n', ...
        summaryTable.method{i}, num_to_str_local(summaryTable.Nc(i)), ...
        num_to_str_local(summaryTable.pdeg(i)), num_to_str_local(summaryTable.refine(i)), ...
        num_to_str_local(summaryTable.h(i)), summaryTable.lambda1(i), ...
        summaryTable.abs_err_lambda1(i), summaryTable.l2_error_u1(i));
end
fclose(fid);
end

function cleanup_legacy_outputs_local(outDir)
%Compute legacy outputs.
patterns = {
    'figure*_example2_*.pdf'
    'figure*_example2_*.png'
    'reference_eigenfunction_*.pdf'
    'reference_eigenfunction_*.png'
    };
for i = 1:numel(patterns)
    files = dir(fullfile(outDir, patterns{i}));
    for k = 1:numel(files)
        delete(fullfile(files(k).folder, files(k).name));
    end
end
oldFields = fullfile(outDir, 'fields', 'example2_compare_fields.mat');
if exist(oldFields, 'file')
    delete(oldFields);
end
oldFields = fullfile(outDir, 'fields', 'example2_compare_u1_fields.mat');
if exist(oldFields, 'file')
    delete(oldFields);
end
end

function cmap = example1_reference_colormap_local(n)
%Compute reference colormap.
anchors255 = [
    55, 105, 105;
    140, 190, 170;
    248, 242, 232;
    236, 170, 145;
    196,  85,  60
    ];
cmap = colormap_from_anchors_local(anchors255, n);
end

function cmap = colormap_from_anchors_local(anchors255, n)
%Build a colormap from anchor colors.
A = anchors255 / 255;
k = size(A, 1);
x = linspace(0, 1, k);
xi = linspace(0, 1, n);
cmap = zeros(n, 3);
for j = 1:3
    cmap(:, j) = interp1(x, A(:, j), xi, 'pchip');
end
cmap = min(max(cmap, 0), 1);
end

function s = num_to_str_local(x)
%Compute to str.
if abs(x - round(x)) < 1e-12
    s = sprintf('%d', round(x));
else
    s = sprintf('%.6g', x);
end
end

function ensure_dir_local(pathName)
%Create a directory when it is missing.
if ~exist(pathName, 'dir')
    mkdir(pathName);
end
end

function tf = is_hybrid_run_usable_local(run, solve_mode, Nc, refine, pdeg)
%Check whether a hybrid run has field data.
tf = isstruct(run) && isfield(run, 'meta') && isstruct(run.meta) ...
    && isfield(run.meta, 'solve_mode') && strcmpi(string(run.meta.solve_mode), string(solve_mode)) ...
    && isfield(run.meta, 'Nc') && isequal(double(run.meta.Nc), double(Nc)) ...
    && isfield(run.meta, 'Refinement') && isequal(double(run.meta.Refinement), double(refine)) ...
    && isfield(run.meta, 'pu') && isequal(double(run.meta.pu), double(pdeg)) ...
    && isfield(run, 'uh') && ~isempty(run.uh) ...
    && isfield(run, 'lambda') && numel(run.lambda) >= 1;
end

function tf = is_single_run_usable_local(run, mode, solve_mode, discValue, pdeg, periodicFlag)
%Compute single run usable.
tf = isstruct(run) && isfield(run, 'meta') && isstruct(run.meta) ...
    && isfield(run.meta, 'mode') && strcmpi(string(run.meta.mode), string(mode)) ...
    && isfield(run.meta, 'solve_mode') && strcmpi(string(run.meta.solve_mode), string(solve_mode)) ...
    && isfield(run, 'uh') && ~isempty(run.uh) ...
    && isfield(run, 'lambda') && numel(run.lambda) >= 1;

if ~tf
    return;
end

switch lower(string(mode))
    case "pw"
        tf = tf && isfield(run.meta, 'Nc') && isequal(double(run.meta.Nc), double(discValue));
    case "iga"
        tf = tf ...
            && isfield(run.meta, 'refine') && isequal(double(run.meta.refine), double(discValue)) ...
            && isfield(run.meta, 'pdeg') && isequal(double(run.meta.pdeg), double(pdeg)) ...
            && isfield(run.meta, 'periodic_bc') && logical(run.meta.periodic_bc) == logical(periodicFlag);
    otherwise
        tf = false;
end
end
