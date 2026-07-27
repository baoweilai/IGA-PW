function faceData = build_face_data_3D(nurbs_refine, a, L, N, method)
% Build trace, normal, and coupling data for all six cube faces.

arguments
    nurbs_refine
    a
    L
    N
    method = 'tensor'
end

% Extract the spline data and plane-wave frequencies.
faceNames = {'x-','x+','y-','y+','z-','z+'};
faceData = cell(numel(faceNames),1);

m  = nurbs_refine.m;
n  = nurbs_refine.n;
l  = nurbs_refine.l;
pu = nurbs_refine.pu;
pv = nurbs_refine.pv;
pw = nurbs_refine.pw;

Ubar = nurbs_refine.Ubar;
Vbar = nurbs_refine.Vbar;
Wbar = nurbs_refine.Wbar;

UBreaks = nurbs_refine.UBreaks;
VBreaks = nurbs_refine.VBreaks;
WBreaks = nurbs_refine.WBreaks;

alpha = 2*pi/L;
kvals = -N:N;

% Select the requested face-assembly formulation.
if ~strcmpi(method, 'legacy')
    faceData = build_face_data_tensor_fast_local( ...
        m, n, l, pu, pv, pw, Ubar, Vbar, Wbar, UBreaks, VBreaks, WBreaks, ...
        a, alpha, kvals, faceNames);
    return;
end

% Assemble trace and normal-derivative matrices on all six faces.
for iface = 1:numel(faceNames)
    fname = faceNames{iface};

    switch fname
        case {'x-','x+'}
            fixedParam = double(strcmp(fname,'x+'));   % u = 0 or 1
            fixedCoord = -a + 2*a*fixedParam;
            normalSign = -1 + 2*double(strcmp(fname,'x+'));

            % Build tangential quadrature for an x-normal face.
            [nodes1, w1, span1] = build_1d_face_nodes(VBreaks, Vbar, pv, a);
            [nodes2, w2, span2] = build_1d_face_nodes(WBreaks, Wbar, pw, a);

            nrow = numel(nodes1) * numel(nodes2);
            TI = spalloc(nrow, m*n*l, nrow * 2*(pv+1)*(pw+1));
            GI = spalloc(nrow, m*n*l, nrow * 2*(pv+1)*(pw+1));

            row = 0;
            for j2 = 1:numel(nodes2)
                w = (nodes2(j2) + a)/(2*a);
                Wders = bspbasisDers(Wbar, pw, w, 1);
                Nw  = Wders(1,:)';
                kList = span2(j2)-pw : span2(j2);

                for j1 = 1:numel(nodes1)
                    row = row + 1;

                    v = (nodes1(j1) + a)/(2*a);
                    Vders = bspbasisDers(Vbar, pv, v, 1);
                    Nv  = Vders(1,:)';
                    jList = span1(j1)-pv : span1(j1);

                    Uders = bspbasisDers(Ubar, pu, fixedParam, 1);
                    if fixedParam == 0
                        Nu  = Uders(1,1:2)';
                        DNu = Uders(2,1:2)';
                        iList = [1,2];
                    else
                        Nu  = Uders(1,end-1:end)';
                        DNu = Uders(2,end-1:end)';
                        iList = [m-1,m];
                    end

                    for kk = 1:(pw+1)
                        for jj = 1:(pv+1)
                            for ii = 1:2
                                col = id3(iList(ii), jList(jj), kList(kk), m, n);
                                val = Nu(ii) * Nv(jj) * Nw(kk);

                                % full normal derivative on x-face
                                dnor = normalSign * (DNu(ii) * Nv(jj) * Nw(kk)) / (2*a);

                                TI(row,col) = TI(row,col) + val;
                                GI(row,col) = GI(row,col) + dnor;
                            end
                        end
                    end
                end
            end

            % Precompute tangential plane-wave phase matrices.
            E1 = exp(1i * alpha * (nodes1(:) * kvals));
            E2 = exp(1i * alpha * (nodes2(:) * kvals));

            % Store the x-face quadrature and coupling operators.
            F = struct();
            F.name       = fname;
            F.type       = 'x';
            F.fixedCoord = fixedCoord;
            F.normalSign = normalSign;
            F.TI         = TI;
            F.GI         = GI;
            F.w          = kron(w2(:), w1(:));
            F.E1         = E1;
            F.E2         = E2;
            F.nq1        = numel(nodes1);
            F.nq2        = numel(nodes2);
            faceData{iface} = F;

        case {'y-','y+'}
            fixedParam = double(strcmp(fname,'y+'));   % v = 0 or 1
            fixedCoord = -a + 2*a*fixedParam;
            normalSign = -1 + 2*double(strcmp(fname,'y+'));

            % Build tangential quadrature for a y-normal face.
            [nodes1, w1, span1] = build_1d_face_nodes(UBreaks, Ubar, pu, a);
            [nodes2, w2, span2] = build_1d_face_nodes(WBreaks, Wbar, pw, a);

            nrow = numel(nodes1) * numel(nodes2);
            TI = spalloc(nrow, m*n*l, nrow * 2*(pu+1)*(pw+1));
            GI = spalloc(nrow, m*n*l, nrow * 2*(pu+1)*(pw+1));

            row = 0;
            for j2 = 1:numel(nodes2)
                w = (nodes2(j2) + a)/(2*a);
                Wders = bspbasisDers(Wbar, pw, w, 1);
                Nw  = Wders(1,:)';
                kList = span2(j2)-pw : span2(j2);

                for j1 = 1:numel(nodes1)
                    row = row + 1;

                    u = (nodes1(j1) + a)/(2*a);
                    Uders = bspbasisDers(Ubar, pu, u, 1);
                    Nu  = Uders(1,:)';
                    iList = span1(j1)-pu : span1(j1);

                    Vders = bspbasisDers(Vbar, pv, fixedParam, 1);
                    if fixedParam == 0
                        Nv  = Vders(1,1:2)';
                        DNv = Vders(2,1:2)';
                        jList = [1,2];
                    else
                        Nv  = Vders(1,end-1:end)';
                        DNv = Vders(2,end-1:end)';
                        jList = [n-1,n];
                    end

                    for kk = 1:(pw+1)
                        for ii = 1:(pu+1)
                            for jj = 1:2
                                col = id3(iList(ii), jList(jj), kList(kk), m, n);
                                val = Nu(ii) * Nv(jj) * Nw(kk);

                                % full normal derivative on y-face
                                dnor = normalSign * (Nu(ii) * DNv(jj) * Nw(kk)) / (2*a);

                                TI(row,col) = TI(row,col) + val;
                                GI(row,col) = GI(row,col) + dnor;
                            end
                        end
                    end
                end
            end

            % Precompute tangential plane-wave phase matrices.
            E1 = exp(1i * alpha * (nodes1(:) * kvals));
            E2 = exp(1i * alpha * (nodes2(:) * kvals));

            % Store the y-face quadrature and coupling operators.
            F = struct();
            F.name       = fname;
            F.type       = 'y';
            F.fixedCoord = fixedCoord;
            F.normalSign = normalSign;
            F.TI         = TI;
            F.GI         = GI;
            F.w          = kron(w2(:), w1(:));
            F.E1         = E1;
            F.E2         = E2;
            F.nq1        = numel(nodes1);
            F.nq2        = numel(nodes2);
            faceData{iface} = F;

        case {'z-','z+'}
            fixedParam = double(strcmp(fname,'z+'));   % w = 0 or 1
            fixedCoord = -a + 2*a*fixedParam;
            normalSign = -1 + 2*double(strcmp(fname,'z+'));

            % Build tangential quadrature for a z-normal face.
            [nodes1, w1, span1] = build_1d_face_nodes(UBreaks, Ubar, pu, a);
            [nodes2, w2, span2] = build_1d_face_nodes(VBreaks, Vbar, pv, a);

            nrow = numel(nodes1) * numel(nodes2);
            TI = spalloc(nrow, m*n*l, nrow * 2*(pu+1)*(pv+1));
            GI = spalloc(nrow, m*n*l, nrow * 2*(pu+1)*(pv+1));

            row = 0;
            for j2 = 1:numel(nodes2)
                v = (nodes2(j2) + a)/(2*a);
                Vders = bspbasisDers(Vbar, pv, v, 1);
                Nv  = Vders(1,:)';
                jList = span2(j2)-pv : span2(j2);

                for j1 = 1:numel(nodes1)
                    row = row + 1;

                    u = (nodes1(j1) + a)/(2*a);
                    Uders = bspbasisDers(Ubar, pu, u, 1);
                    Nu  = Uders(1,:)';
                    iList = span1(j1)-pu : span1(j1);

                    Wders = bspbasisDers(Wbar, pw, fixedParam, 1);
                    if fixedParam == 0
                        Nw  = Wders(1,1:2)';
                        DNw = Wders(2,1:2)';
                        kList = [1,2];
                    else
                        Nw  = Wders(1,end-1:end)';
                        DNw = Wders(2,end-1:end)';
                        kList = [l-1,l];
                    end

                    for jj = 1:(pv+1)
                        for ii = 1:(pu+1)
                            for kk = 1:2
                                col = id3(iList(ii), jList(jj), kList(kk), m, n);
                                val = Nu(ii) * Nv(jj) * Nw(kk);

                                % full normal derivative on z-face
                                dnor = normalSign * (Nu(ii) * Nv(jj) * DNw(kk)) / (2*a);

                                TI(row,col) = TI(row,col) + val;
                                GI(row,col) = GI(row,col) + dnor;
                            end
                        end
                    end
                end
            end

            % Precompute tangential plane-wave phase matrices.
            E1 = exp(1i * alpha * (nodes1(:) * kvals));
            E2 = exp(1i * alpha * (nodes2(:) * kvals));

            % Store the z-face quadrature and coupling operators.
            F = struct();
            F.name       = fname;
            F.type       = 'z';
            F.fixedCoord = fixedCoord;
            F.normalSign = normalSign;
            F.TI         = TI;
            F.GI         = GI;
            F.w          = kron(w2(:), w1(:));
            F.E1         = E1;
            F.E2         = E2;
            F.nq1        = numel(nodes1);
            F.nq2        = numel(nodes2);
            faceData{iface} = F;

        otherwise
            error('Unknown face name.');
    end
end

end

function [nodes, wphys, span] = build_1d_face_nodes(Breaks, Ubar, p, a)
% Build one-dimensional face quadrature nodes, weights, and spans.
[gp, gw] = grule(10*p + 5);
nel = length(Breaks) - 1;
nq  = numel(gp);

nodes = zeros(nel*nq, 1);
wphys = zeros(nel*nq, 1);
span  = zeros(nel*nq, 1);

% Map each knot-span quadrature rule to physical coordinates.
cnt = 0;
for e = 1:nel
    aa = Breaks(e);
    bb = Breaks(e+1);

    Jp = (bb-aa)/2;
    sp = findspan(Ubar, p, Breaks(e));

    for i = 1:nq
        cnt = cnt + 1;
        xi = ((bb-aa)*gp(i) + aa + bb)/2;      % param coordinate
        nodes(cnt) = -a + 2*a*xi;              % physical coordinate on tangent line
        wphys(cnt) = Jp * (2*a) * gw(i);       % physical line weight
        span(cnt)  = sp;
    end
end
end

function idx = id3(i,j,k,m,n)
% Map three-dimensional control indices to one vector index.
idx = i + (j-1)*m + (k-1)*m*n;
end

function faceData = build_face_data_tensor_fast_local( ...
m, n, l, pu, pv, pw, Ubar, Vbar, Wbar, UBreaks, VBreaks, WBreaks, ...
    a, alpha, kvals, faceNames)
% Build all face operators by tensor products.
faceData = cell(numel(faceNames), 1);

% Build one-dimensional quadrature basis matrices.
[uNodes, uWeights, Bu] = build_1d_face_basis_matrix_local(UBreaks, Ubar, pu, m, a);
[vNodes, vWeights, Bv] = build_1d_face_basis_matrix_local(VBreaks, Vbar, pv, n, a);
[wNodes, wWeights, Bw] = build_1d_face_basis_matrix_local(WBreaks, Wbar, pw, l, a);

% Form Kronecker trace and derivative operators for each face.
for iface = 1:numel(faceNames)
    fname = faceNames{iface};
    fixedParam = double(fname(2) == '+');
    fixedCoord = -a + 2 * a * fixedParam;
    normalSign = -1 + 2 * double(fname(2) == '+');

    switch fname(1)
        case 'x'
            [Bfix, Dfix] = build_endpoint_basis_local(Ubar, pu, m, fixedParam, normalSign, a);
            TI = kron(Bw, kron(Bv, Bfix));
            GI = kron(Bw, kron(Bv, Dfix));
            nodes1 = vNodes;
            nodes2 = wNodes;
            weights1 = vWeights;
            weights2 = wWeights;
            B1 = Bv;
            B2 = Bw;
            ftype = 'x';

        case 'y'
            [Bfix, Dfix] = build_endpoint_basis_local(Vbar, pv, n, fixedParam, normalSign, a);
            TI = kron(Bw, kron(Bfix, Bu));
            GI = kron(Bw, kron(Dfix, Bu));
            nodes1 = uNodes;
            nodes2 = wNodes;
            weights1 = uWeights;
            weights2 = wWeights;
            B1 = Bu;
            B2 = Bw;
            ftype = 'y';

        case 'z'
            [Bfix, Dfix] = build_endpoint_basis_local(Wbar, pw, l, fixedParam, normalSign, a);
            TI = kron(Bfix, kron(Bv, Bu));
            GI = kron(Dfix, kron(Bv, Bu));
            nodes1 = uNodes;
            nodes2 = vNodes;
            weights1 = uWeights;
            weights2 = vWeights;
            B1 = Bu;
            B2 = Bv;
            ftype = 'z';

        otherwise
            error('Unknown face name.');
    end

    % Store the tensor factors and plane-wave phases.
    F = struct();
    F.name = fname;
    F.type = ftype;
    F.fixedCoord = fixedCoord;
    F.normalSign = normalSign;
    F.TI = TI;
    F.GI = GI;
    F.w = kron(weights2(:), weights1(:));
    F.E1 = exp(1i * alpha * (nodes1(:) * kvals));
    F.E2 = exp(1i * alpha * (nodes2(:) * kvals));
    F.nq1 = numel(nodes1);
    F.nq2 = numel(nodes2);
    F.assembly_method = 'tensor_fast';
    F.tensor = struct( ...
        'B1', B1, ...
        'B2', B2, ...
        'Bfix', Bfix, ...
        'Dfix', Dfix, ...
        'w1', weights1(:), ...
        'w2', weights2(:));
    faceData{iface} = F;
end
end

function [nodes, weights, B] = build_1d_face_basis_matrix_local(Breaks, Ubar, p, nCtrl, a)
% Build a one-dimensional face quadrature basis matrix.
[gp, gw] = grule(10 * p + 5);
nel = numel(Breaks) - 1;
nq = numel(gp);
nRows = nel * nq;
nEntries = nRows * (p + 1);

nodes = zeros(nRows, 1);
weights = zeros(nRows, 1);
rows = zeros(nEntries, 1);
cols = zeros(nEntries, 1);
vals = zeros(nEntries, 1);

% Fill the sparse basis rows span by span.
row = 0;
cursor = 0;
for e = 1:nel
    aa = Breaks(e);
    bb = Breaks(e + 1);
    Jp = (bb - aa) / 2;
    span = findspan(Ubar, p, Breaks(e));
    active = span - p : span;

    for q = 1:nq
        row = row + 1;
        xi = ((bb - aa) * gp(q) + aa + bb) / 2;
        ders = bspbasisDers(Ubar, p, xi, 1);

        nodes(row) = -a + 2 * a * xi;
        weights(row) = Jp * (2 * a) * gw(q);

        loc = cursor + (1:(p + 1));
        rows(loc) = row;
        cols(loc) = active(:);
        vals(loc) = ders(1, :).';
        cursor = cursor + p + 1;
    end
end

B = sparse(rows, cols, vals, nRows, nCtrl);
end

function [Bfix, Dfix] = build_endpoint_basis_local(Ubar, p, nCtrl, fixedParam, normalSign, a)
% Build endpoint trace and outward-normal derivative rows.
ders = bspbasisDers(Ubar, p, fixedParam, 1);
if fixedParam == 0
    idx = 1:min(2, nCtrl);
    src = 1:numel(idx);
else
    idx = max(1, nCtrl - 1):nCtrl;
    src = size(ders, 2) - numel(idx) + 1:size(ders, 2);
end

Bfix = sparse(1, idx, ders(1, src), 1, nCtrl);
Dfix = sparse(1, idx, normalSign * ders(2, src) / (2 * a), 1, nCtrl);
end
