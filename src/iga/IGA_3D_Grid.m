function nurbsInfo = IGA_3D_Grid(knotU, knotV, knotW, pu, pv, pw, Refinement)
% Build the 3-D IGA mesh data.

% Refine the knot vectors and store the global basis sizes.
[Ubar, Vbar, Wbar, dof] = IGAknotRefineVolume(knotU, pu, knotV, pv, knotW, pw, Refinement);

nurbsInfo.Ubar = Ubar;
nurbsInfo.Vbar = Vbar;
nurbsInfo.Wbar = Wbar;
nurbsInfo.dof  = dof;

UBreaks = unique(Ubar);
VBreaks = unique(Vbar);
WBreaks = unique(Wbar);

nurbsInfo.UBreaks = UBreaks;
nurbsInfo.VBreaks = VBreaks;
nurbsInfo.WBreaks = WBreaks;

m = length(Ubar) - pu - 1;
n = length(Vbar) - pv - 1;
l = length(Wbar) - pw - 1;

nurbsInfo.m = m;
nurbsInfo.n = n;
nurbsInfo.l = l;

nurbsInfo.n_dofs         = m*n*l;
nurbsInfo.n_dofs_domains = m*n*l;

nurbsInfo.pu = pu;
nurbsInfo.pv = pv;
nurbsInfo.pw = pw;

% Build the element intervals and counts in each direction.
uNoEs = length(UBreaks) - 1;
vNoEs = length(VBreaks) - 1;
wNoEs = length(WBreaks) - 1;
NoEs  = uNoEs * vNoEs * wNoEs;

nurbsInfo.uNoEs = uNoEs;
nurbsInfo.vNoEs = vNoEs;
nurbsInfo.wNoEs = wNoEs;
nurbsInfo.NoEs  = NoEs;

n_ele_dofs = (pu+1)*(pv+1)*(pw+1);

% Build element connectivity and parametric bounds.
Element       = zeros(NoEs, n_ele_dofs);
Coordinate    = zeros(NoEs, 6);   % [u1 u2 v1 v2 w1 w2]
knotSpanIndex = zeros(NoEs, 3);   % [ispan jspan kspan]

e = 0;
for k1 = 1:wNoEs
    for j1 = 1:vNoEs
        for i1 = 1:uNoEs
            e = e + 1;

            ue = [UBreaks(i1), UBreaks(i1+1)];
            ve = [VBreaks(j1), VBreaks(j1+1)];
            we = [WBreaks(k1), WBreaks(k1+1)];

            Coordinate(e,:) = [ue, ve, we];

            ispan = findspan(Ubar, pu, UBreaks(i1));
            jspan = findspan(Vbar, pv, VBreaks(j1));
            kspan = findspan(Wbar, pw, WBreaks(k1));

            knotSpanIndex(e,:) = [ispan, jspan, kspan];

            row = zeros(1, n_ele_dofs);
            cnt = 0;

            % Order local basis indices with u fastest, then v and w.
            for kk = 0:pw
                for jj = 0:pv
                    for ii = 0:pu
                        cnt = cnt + 1;

                        gi = ispan - pu + ii;
                        gj = jspan - pv + jj;
                        gk = kspan - pw + kk;

                        row(cnt) = gi + (gj-1)*m + (gk-1)*m*n;
                    end
                end
            end

            Element(e,:) = row;
        end
    end
end

% Store the element tables in the mesh structure.
nurbsInfo.Element       = Element;
nurbsInfo.Coordinate    = Coordinate;
nurbsInfo.knotSpanIndex = knotSpanIndex;

end
