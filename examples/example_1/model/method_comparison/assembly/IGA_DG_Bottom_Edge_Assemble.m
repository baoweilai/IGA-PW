function [P, S] = IGA_DG_Bottom_Edge_Assemble( ...
nurbs_original, nurbs_refine, pw_index, plane_wave_dofs_index, L, n_dofs)
% Assemble DG interface matrices on the bottom edge:
%   P : jump-jump term
%   S : jump-average term

DIM = 2;

%% Input NURBS patch data
ConPts_o  = nurbs_original.ConPts;
weights_o = nurbs_original.weights;
knotU_o   = nurbs_original.knotU;
knotV_o   = nurbs_original.knotV;
pu_o      = nurbs_original.pu;
pv_o      = nurbs_original.pv;

%% Refined NURBS data
knotU = nurbs_refine.Ubar;
knotV = nurbs_refine.Vbar;
pu    = nurbs_refine.pu;
pv    = nurbs_refine.pv;

%% Initialize global matrices
S = sparse(n_dofs, n_dofs);
P = S;

% Plus side normal points from inner domain to outer domain
UBreaks          = nurbs_refine.UBreaks;
uNoEs            = length(UBreaks) - 1;
bottom_edge_dofs = nurbs_refine.bottom_edge_dofs;

%% Parametric intervals of bottom edges
bottom_edge_node = zeros(uNoEs, 2);
for e = 1:uNoEs
    bottom_edge_node(e, :) = [UBreaks(e), UBreaks(e + 1)];
end

%% Gaussian quadrature on [-1,1]
[gp, gw] = grule(10 * pu + 5);
n_gp = length(gp);

%% Bottom interface: v = 0
v_bottom_inner = 0;

%% Plane wave data
edge_dofs_minus  = plane_wave_dofs_index;
n_pw_basis       = size(pw_index, 1);
basis_grad_minus = zeros(n_pw_basis, DIM);
basis_minus      = zeros(n_pw_basis, 1);
Omega_area       = L * L;

%% Loop over bottom edges
for e = 1:uNoEs
    edge_dofs_plus = bottom_edge_dofs(e, :);
    edge_dofs      = [edge_dofs_plus, edge_dofs_minus];
    n_edge_dofs    = length(edge_dofs);

    ue = bottom_edge_node(e, :);
    a  = ue(1);
    b  = ue(2);
    J1 = (b - a) / 2;  % Map [-1,1] to [a,b]

    edge_jump_Ae    = zeros(n_edge_dofs, n_edge_dofs);
    edge_average_Ae = zeros(n_edge_dofs, n_edge_dofs);

    %% Quadrature on current edge
    for i = 1:n_gp
        u = ((b - a) * gp(i) + a + b) / 2;

% Geometry and Jacobian on the input patch
        [F, DF_plus] = NurbsSurface( ...
            ConPts_o, weights_o, knotU_o, pu_o, u, knotV_o, pv_o, v_bottom_inner);

        tau    = DF_plus(:, 1);          % Tangent vector
        ds     = norm(tau);              % Arc-length factor
        normal = [tau(2); -tau(1)] / ds; % Unit normal
        Jacobi = J1 * gw(i) * ds;        % Edge quadrature weight

        %% Plus-side NURBS basis
        Uders_plus = bspbasisDers(knotU, pu, u, 1);
        Nu_plus    = Uders_plus(1, :)';
        DNu_plus   = Uders_plus(2, :)';

        Vders_plus = bspbasisDers(knotV, pv, v_bottom_inner, 1);
        Nv_plus    = Vders_plus(1, 1:2);
        DNv_plus   = Vders_plus(2, 1:2);

        basis_plus = Nu_plus * Nv_plus;
        basis_plus = basis_plus(:);

        DNu_v_plus = DNu_plus * Nv_plus;
        DNu_v_plus = DNu_v_plus(:);

        DNv_u_plus = Nu_plus * DNv_plus;
        DNv_u_plus = DNv_u_plus(:);

        % Physical gradients of NURBS basis
        basis_grad_plus = [DNu_v_plus, DNv_u_plus] / DF_plus;

        %% Minus-side plane wave basis
        for k = 1:n_pw_basis
            basis_minus(k) = exp(-1i * 2 * pi / L * pw_index(k, :) * F) ...
                / sqrt(Omega_area);

            basis_grad_minus(k, :) = (1i * 2 * pi / L) ...
                * pw_index(k, :) * basis_minus(k)';
        end

        %% Jump term
        edge_jump = [basis_plus; -basis_minus];
        edge_jump_Ae = edge_jump_Ae + edge_jump * (edge_jump') * Jacobi;

        %% Average term
        edge_average = [basis_grad_plus; basis_grad_minus] * normal / 2;
        edge_average = edge_average.';  % transpose only
        edge_average_Ae = edge_average_Ae + edge_jump * edge_average * Jacobi;
    end

    %% Assemble into global matrices
    P(edge_dofs, edge_dofs) = P(edge_dofs, edge_dofs) + edge_jump_Ae;
    S(edge_dofs, edge_dofs) = S(edge_dofs, edge_dofs) + edge_average_Ae;
end

end
