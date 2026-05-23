function [Ubar,Vbar,Wbar,dof] = IGAknotRefineVolume(knotU,pu,knotV,pv,knotW,pw,Refinement)
%Refine NURBS volume knot vectors.

if isnumeric(Refinement) && isscalar(Refinement)
    Ubar = knotU;
    Vbar = knotV;
    Wbar = knotW;

    for i = 1:Refinement
        UBreks = unique(Ubar);
        VBreks = unique(Vbar);
        WBreks = unique(Wbar);

        Xu = (UBreks(1:end-1) + UBreks(2:end)) / 2;
        Xv = (VBreks(1:end-1) + VBreks(2:end)) / 2;
        Xw = (WBreks(1:end-1) + WBreks(2:end)) / 2;

        Ubar = sort([Ubar, Xu]);
        Vbar = sort([Vbar, Xv]);
        Wbar = sort([Wbar, Xw]);
    end
elseif isstruct(Refinement) && isfield(Refinement, 'mode') && strcmpi(Refinement.mode, 'nelem')
    nElem = Refinement.value;

    if isscalar(nElem)
        nElemU = nElem;
        nElemV = nElem;
        nElemW = nElem;
    else
        assert(numel(nElem) == 3, ...
            'Refinement.value must be scalar or 3-vector in nelem mode.');
        nElemU = nElem(1);
        nElemV = nElem(2);
        nElemW = nElem(3);
    end

    assert(nElemU >= 1 && nElemV >= 1 && nElemW >= 1, ...
        'Each number of elements must be >= 1.');

    Ubar = build_open_uniform_knot(knotU(1), knotU(end), pu, nElemU);
    Vbar = build_open_uniform_knot(knotV(1), knotV(end), pv, nElemV);
    Wbar = build_open_uniform_knot(knotW(1), knotW(end), pw, nElemW);
else
    error(['Unsupported Refinement input.' newline ...
        'Use either a scalar refinement count or ' ...
        'struct(''mode'',''nelem'',''value'',nElem).']);
end

nu = length(Ubar) - pu - 1;
nv = length(Vbar) - pv - 1;
nw = length(Wbar) - pw - 1;
dof = nu * nv * nw;
end

function K = build_open_uniform_knot(a, b, p, nElem)
%Build an open uniform knot vector on [a,b].
if nElem == 1
    K = [repmat(a, 1, p + 1), repmat(b, 1, p + 1)];
else
    breaks = linspace(a, b, nElem + 1);
    interior = breaks(2:end-1);
    K = [repmat(a, 1, p + 1), interior, repmat(b, 1, p + 1)];
end
end
