function [supported, support] = build_affine_cube_fast_support_3D( ...
nurbs_original, nurbs_refine, n_gp)
% Precompute tensor-product data for affine-cube fast assembly paths.

support = struct();
[supported, geom] = detect_affine_cube_geometry_local(nurbs_original);
if ~supported
    return;
end

[gp, gw] = grule(n_gp);

support.geom = geom;
support.gp = gp;
support.gw = gw;
support.Element = nurbs_refine.Element;
support.Coordinate = nurbs_refine.Coordinate;
support.NoEs = nurbs_refine.NoEs;
support.n_dofs = nurbs_refine.n_dofs;
support.uNoEs = nurbs_refine.uNoEs;
support.vNoEs = nurbs_refine.vNoEs;
support.wNoEs = nurbs_refine.wNoEs;
support.pu = nurbs_refine.pu;
support.pv = nurbs_refine.pv;
support.pw = nurbs_refine.pw;
support.n_ele_dofs = (support.pu + 1) * (support.pv + 1) * (support.pw + 1);

support.u_cache = precompute_basis_line_local( ...
    nurbs_refine.UBreaks, nurbs_refine.Ubar, support.pu, gp, gw, geom.x0, geom.sx);
support.v_cache = precompute_basis_line_local( ...
    nurbs_refine.VBreaks, nurbs_refine.Vbar, support.pv, gp, gw, geom.y0, geom.sy);
support.w_cache = precompute_basis_line_local( ...
    nurbs_refine.WBreaks, nurbs_refine.Wbar, support.pw, gp, gw, geom.z0, geom.sz);
end

function cache = precompute_basis_line_local(Breaks, knot, p, gp, gw, phys0, physScale)
%Compute basis line.
n_ele = numel(Breaks) - 1;
nq = numel(gp);
cache = cell(n_ele, 1);

for e = 1:n_ele
    interval = Breaks(e:e+1);
    J = (interval(2) - interval(1)) / 2;
    param_pts = ((interval(2) - interval(1)) * gp + interval(1) + interval(2)) / 2;

    N = zeros(p + 1, nq);
    D_phys = zeros(p + 1, nq);
    for iq = 1:nq
        ders = bspbasisDers(knot, p, param_pts(iq), 1);
        N(:, iq) = ders(1, :)';
        D_phys(:, iq) = ders(2, :)' / physScale;
    end

    entry = struct();
    entry.N = N;
    entry.D_phys = D_phys;
    entry.w_line = J * gw(:);
    entry.phys_pts = phys0 + physScale * param_pts(:);
    cache{e} = entry;
end
end

function [tf, geom] = detect_affine_cube_geometry_local(nurbs_original)
%Detect affine geometry for a cube patch.
geom = struct('x0', 0, 'y0', 0, 'z0', 0, 'sx', 0, 'sy', 0, 'sz', 0, 'abs_detDF', 0);
tf = false;

ConPts = nurbs_original.ConPts;
weights = nurbs_original.weights;
[nu, nv, nw, dim] = size(ConPts);
if dim ~= 3 || nu < 2 || nv < 2 || nw < 2
    return;
end
if any(abs(weights(:) - 1) > 1e-14)
    return;
end

x_grid = ConPts(:, :, :, 1);
y_grid = ConPts(:, :, :, 2);
z_grid = ConPts(:, :, :, 3);

if max_abs_diff_local(x_grid, 2) > 1e-12 || max_abs_diff_local(x_grid, 3) > 1e-12
    return;
end
if max_abs_diff_local(y_grid, 1) > 1e-12 || max_abs_diff_local(y_grid, 3) > 1e-12
    return;
end
if max_abs_diff_local(z_grid, 1) > 1e-12 || max_abs_diff_local(z_grid, 2) > 1e-12
    return;
end

sx = x_grid(end, 1, 1) - x_grid(1, 1, 1);
sy = y_grid(1, end, 1) - y_grid(1, 1, 1);
sz = z_grid(1, 1, end) - z_grid(1, 1, 1);
if abs(sx) <= 1e-14 || abs(sy) <= 1e-14 || abs(sz) <= 1e-14
    return;
end

geom.x0 = x_grid(1, 1, 1);
geom.y0 = y_grid(1, 1, 1);
geom.z0 = z_grid(1, 1, 1);
geom.sx = sx;
geom.sy = sy;
geom.sz = sz;
geom.abs_detDF = abs(sx * sy * sz);

tf = true;
end

function v = max_abs_diff_local(A, dim)
%Compute the maximum absolute difference.
tmp = abs(diff(A, 1, dim));
if isempty(tmp)
    v = 0;
else
    v = max(tmp(:));
end
end
