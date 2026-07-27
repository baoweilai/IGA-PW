function [Pfun, info] = tbprec_rfp( ...
    PGamma, pwIdx, diagShifted, nTotal, traceFactorBuilder, cfg, caseInfo)

% Build the RFP trace preconditioner.
assert(cfg.patternTol == 1e-12, 'TB-DG patternTol must equal 1e-12.');
assert(cfg.targetShift == 0, 'TB-DG targetShift must equal zero.');
assert(exist('rfp_chol_mex', 'file') == 3, ...
    'Compile rfp_chol_mex before building the RFP preconditioner.');

patternTol = cfg.patternTol;
targetShift = cfg.targetShift;

Ppat = PGamma;
[ii, jj, vv] = find(Ppat);
keep = abs(vv) >= patternTol;
Ppat = sparse(ii(keep), jj(keep), vv(keep), nTotal, nTotal);
[ii, jj] = find(Ppat);

gamma = unique([pwIdx(:); ii; jj]);
eta = setdiff((1:nTotal).', gamma);

d = abs(diagShifted(eta));
if any(~isfinite(d)) || any(d == 0)
    error(['TB-DG interior diagonal failed: K=%g, p=%g, nelem=%g, ' ...
        'shift=%.17g.'], caseInfo.K, caseInfo.p, caseInfo.nelem, targetShift);
end

rfp_chol_mex('free');
traceInfo = traceFactorBuilder(gamma);
Pfun = @(x) apply_tbprec_rfp_local(x, gamma, eta, d);

info = traceInfo;
info.patternTol = patternTol;
info.targetShift = targetShift;
info.gamma = gamma;
info.eta = eta;
info.nGamma = numel(gamma);
info.nEta = numel(eta);
info.factorType = 'rfp_zpftrf';
end

function y = apply_tbprec_rfp_local(x, gamma, eta, d)
% Apply the RFP trace preconditioner.
y = zeros(size(x), 'like', x);
y(eta, :) = x(eta, :) ./ d;
y(gamma, :) = rfp_chol_mex('solve', x(gamma, :));
end
