function E = build_iga_grid_eval_matrix_3D(nurbs_refine, X, Y, Z, a)
%Build IGA grid eval matrix 3D.

U = nurbs_refine.Ubar;
V = nurbs_refine.Vbar;
W = nurbs_refine.Wbar;
pu = nurbs_refine.pu;
pv = nurbs_refine.pv;
pw = nurbs_refine.pw;
m = nurbs_refine.m;
n = nurbs_refine.n;
nDofs = nurbs_refine.n_dofs;

nq = numel(X);
nzPerRow = (pu + 1) * (pv + 1) * (pw + 1);
rows = zeros(nq * nzPerRow, 1);
cols = zeros(nq * nzPerRow, 1);
vals = zeros(nq * nzPerRow, 1);

cursor = 0;
for q = 1:nq
    u = min(max((X(q) + a) / (2 * a), 0), 1);
    v = min(max((Y(q) + a) / (2 * a), 0), 1);
    w = min(max((Z(q) + a) / (2 * a), 0), 1);

    ispan = findspan(U, pu, u);
    jspan = findspan(V, pv, v);
    kspan = findspan(W, pw, w);

    Nu = bspbasisDers(U, pu, u, 1);
    Nv = bspbasisDers(V, pv, v, 1);
    Nw = bspbasisDers(W, pw, w, 1);
    Nu = Nu(1, :)';
    Nv = Nv(1, :)';
    Nw = Nw(1, :)';

    ii = ispan - pu:ispan;
    jj = jspan - pv:jspan;
    kk = kspan - pw:kspan;

    for kz = 1:(pw + 1)
        for jy = 1:(pv + 1)
            for ix = 1:(pu + 1)
                cursor = cursor + 1;
                rows(cursor) = q;
                cols(cursor) = ii(ix) + (jj(jy) - 1) * m + (kk(kz) - 1) * m * n;
                vals(cursor) = Nu(ix) * Nv(jy) * Nw(kz);
            end
        end
    end
end

E = sparse(rows, cols, vals, nq, nDofs);
end
