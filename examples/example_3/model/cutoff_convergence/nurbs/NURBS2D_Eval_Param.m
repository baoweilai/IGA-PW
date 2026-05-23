function [ids, R, dRdx] = NURBS2D_Eval_Param(refine, x, y)
%Evaluate a 2-D NURBS parameter map.

U  = refine.Ubar;    V  = refine.Vbar;
p  = refine.pu;      q  = refine.pv;
m  = length(U) - p - 1;
n  = length(V) - q - 1;

xmin = refine.bbox(1,1); xmax = refine.bbox(1,2);
ymin = refine.bbox(2,1); ymax = refine.bbox(2,2);
Lx   = xmax - xmin;       Ly   = ymax - ymin;

u = (x - xmin)/Lx;
v = (y - ymin)/Ly;

u = min(max(u, U(p+1)), U(end-p));
v = min(max(v, V(q+1)), V(end-q));

spanU = findspan(U, p, u);
spanV = findspan(V, q, v);

uders = bspbasisDers(U, p, u, 1);
vders = bspbasisDers(V, q, v, 1);   % 2 x (q+1)

Nu  = uders(1, :).';   DNu = uders(2, :).';
Nv  = vders(1, :);     DNv = vders(2, :);

% R(u,v) = Nu * Nv
R     = (Nu * Nv);    R     = R(:);
DNuNv = (DNu * Nv);   DNuNv = DNuNv(:);
NuDNv = (Nu  * DNv);  NuDNv = NuDNv(:);

dRdx        = zeros(numel(R), 2);
dRdx(:, 1)  = DNuNv / Lx;
dRdx(:, 2)  = NuDNv / Ly;

ids = zeros(numel(R), 1);
k = 0;
for j = (spanV - q) : spanV
    for i = (spanU - p) : spanU
        k = k + 1;
        ids(k) = i + (j - 1) * m;  % 1-based
    end
end
end
