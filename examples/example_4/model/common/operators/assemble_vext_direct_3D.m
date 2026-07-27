function [Vmat, meta] = assemble_vext_direct_3D( ...
    nurbsOriginal, nurbsRefine, kVr, nVr, ewald, nGauss, opts)
% Assemble the Ewald external potential at IGA Gauss points.

arguments
    nurbsOriginal
    nurbsRefine
    kVr
    nVr (1,1) double
    ewald struct
    nGauss (1,1) double
    opts struct
end

% Validate affine-cube support and build reciprocal coefficients.
assert(opts.use_affine_cube_fast, ...
    'Direct Ewald assembly requires the affine-cube fast path.');
[supported, support] = build_affine_cube_fast_support_3D( ...
    nurbsOriginal, nurbsRefine, nGauss);
assert(supported, 'The current IGA patch is not an affine cube.');

tTotal = tic;
[G, coeff] = reciprocal_data_local(kVr, nVr, ewald);

Element = support.Element;
NoEs = support.NoEs;
nDofs = support.n_dofs;
uNoEs = support.uNoEs;
vNoEs = support.vNoEs;
wNoEs = support.wNoEs;
nElementDofs = support.n_ele_dofs;
nLocalEntries = nElementDofs*nElementDofs;

% Allocate sparse triplets and timing accumulators.
Vval = zeros(NoEs*nLocalEntries, 1);
Iind = zeros(NoEs*nLocalEntries, 1);
Jind = zeros(NoEs*nLocalEntries, 1);
gidx = 1;

tEvaluate = 0;
tLocal = 0;
tScatter = 0;

% Evaluate and integrate the potential over each tensor-product element.
for kw = 1:wNoEs
    wdata = support.w_cache{kw};
    for jv = 1:vNoEs
        vdata = support.v_cache{jv};
        for iu = 1:uNoEs
            e = iu+(jv-1)*uNoEs+(kw-1)*uNoEs*vNoEs;
            udata = support.u_cache{iu};
            row = Element(e, :);

            t0 = tic;
            [Phi, wJ, xPts, yPts, zPts] = element_data_local( ...
                udata, vdata, wdata, support.geom.abs_detDF);
            Vq = evaluate_ewald_local(ewald, G, coeff, xPts, yPts, zPts);
            tEvaluate = tEvaluate+toc(t0);

            t0 = tic;
            weightedPhi = bsxfun(@times, Phi, reshape(wJ.*Vq, 1, []));
            Vloc = weightedPhi*Phi';
            tLocal = tLocal+toc(t0);

            t0 = tic;
            idx = gidx:(gidx+nLocalEntries-1);
            rowVector = row(:);
            Iind(idx) = repmat(rowVector, nElementDofs, 1);
            Jind(idx) = repelem(rowVector, nElementDofs, 1);
            Vval(idx) = Vloc(:);
            gidx = gidx+nLocalEntries;
            tScatter = tScatter+toc(t0);
        end
    end
end

% Assemble the global potential matrix and timing data.
Vmat = sparse(Iind, Jind, Vval, nDofs, nDofs);
Vmat = 0.5*(Vmat+Vmat');

meta = struct();
meta.method = 'iga_gauss_direct_ewald';
meta.t_evaluate = tEvaluate;
meta.t_local = tLocal;
meta.t_scatter = tScatter;
meta.t_total = toc(tTotal);
meta.n_reciprocal_modes = numel(coeff);
end

function [G, coeff] = reciprocal_data_local(kVr, nVr, ewald)
% Build reciprocal vectors and Ewald coefficients.
modes = kVr(1:nVr, :);
positiveRepresentative = modes(:, 1) > 0 ...
    | (modes(:, 1) == 0 & modes(:, 2) > 0) ...
    | (modes(:, 1) == 0 & modes(:, 2) == 0 & modes(:, 3) > 0);
G = (2*pi/ewald.L)*modes(positiveRepresentative, :);
g2 = sum(G.^2, 2);
coeff = -2*ewald.charge*(4*pi/ewald.L^3) ...
    .* exp(-g2/(4*ewald.mu^2))./g2;
end

function Vq = evaluate_ewald_local(ewald, G, coeff, xPts, yPts, zPts)
% Evaluate the Ewald potential at quadrature points.
r = sqrt(xPts.^2+yPts.^2+zPts.^2);
r = max(r, 1e-14);
Vq = -ewald.charge*erfc(ewald.mu*r)./r;
if ~isempty(coeff)
    phase = xPts*G(:, 1).'+yPts*G(:, 2).'+zPts*G(:, 3).';
    Vq = Vq+cos(phase)*coeff;
end
Vq = Vq+ewald.charge*2*ewald.mu/sqrt(pi);
end

function [Phi, wJ, xPts, yPts, zPts] = element_data_local( ...
    udata, vdata, wdata, absDetDF)
% Build basis, quadrature, and physical-point data for one element.
Phi = kron(wdata.N, kron(vdata.N, udata.N));
wJ = absDetDF*kron(wdata.w_line, kron(vdata.w_line, udata.w_line));
[X, Y, Z] = ndgrid(udata.phys_pts, vdata.phys_pts, wdata.phys_pts);
xPts = X(:);
yPts = Y(:);
zPts = Z(:);
end
