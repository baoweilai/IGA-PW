function gridCache = build_midpoint_grid_cache_3D(L, mFFT, a)
%Build midpoint grid cache 3D.

dx = L / mFFT;
xmid = -L / 2 + dx / 2 + (0:mFFT-1) * dx;
inner_idx = find(abs(xmid) <= a + 1e-12);
x_inner = xmid(inner_idx);
[Xin, Yin, Zin] = ndgrid(x_inner, x_inner, x_inner);

gridCache = struct();
gridCache.L = L;
gridCache.a = a;
gridCache.dx = dx;
gridCache.mFFT = mFFT;
gridCache.xmid = xmid;
gridCache.inner_idx = inner_idx;
gridCache.n_inner = numel(inner_idx);
gridCache.x_inner = Xin(:);
gridCache.y_inner = Yin(:);
gridCache.z_inner = Zin(:);
end
