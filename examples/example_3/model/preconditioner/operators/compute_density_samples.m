function rhoData = compute_density_samples( ...
u, n_dofs_nurbs, pw_dofs_indices, ...
    nurbs_original, nurbs_refine, n_gp, ...
    k_pw, L, a, m_pw)

% Evaluate the inner and outer density samples.
uI = u(1:n_dofs_nurbs);
uA = u(pw_dofs_indices);

rho_in = sample_rho_nurbs(uI, nurbs_original, nurbs_refine, n_gp);
[rho_pw, dx_pw] = sample_rho_pw_grid(uA, k_pw, L, a, m_pw);

rhoData = struct();
rhoData.rho_in = rho_in;
rhoData.rho_pw = rho_pw;
rhoData.dx_pw  = dx_pw;
end
