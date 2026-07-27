function [P,S] = IGA_DG_Face_Assemble_3D(nurbs_original, nurbs_refine, pw_index, plane_wave_dofs_index, L, n_dofs, faceName)
% Assemble jump and average operators on one three-dimensional interface face.

% Read the NURBS geometry, refined basis, and element breaks.
DIM = 3;

ConPts_o   = nurbs_original.ConPts;
weights_o  = nurbs_original.weights;
knotU_o    = nurbs_original.knotU;
knotV_o    = nurbs_original.knotV;
knotW_o    = nurbs_original.knotW;
pu_o       = nurbs_original.pu;
pv_o       = nurbs_original.pv;
pw_o       = nurbs_original.pw;

knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
knotW      = nurbs_refine.Wbar;
pu         = nurbs_refine.pu;
pv         = nurbs_refine.pv;
pw         = nurbs_refine.pw;

UBreaks = nurbs_refine.UBreaks;
VBreaks = nurbs_refine.VBreaks;
WBreaks = nurbs_refine.WBreaks;

m = nurbs_refine.m;
n = nurbs_refine.n;
l = nurbs_refine.l;

S = sparse(n_dofs, n_dofs);
P = sparse(n_dofs, n_dofs);

% Prepare the plane-wave traces and global interface matrices.
edge_dofs_minus  = plane_wave_dofs_index;
n_pw_basis       = size(pw_index,1);
basis_grad_minus = zeros(n_pw_basis, DIM);
basis_minus      = zeros(n_pw_basis, 1);
Omega_vol        = L^3;

% Determine tangential directions on the selected face
switch faceName
    case {'x-','x+'}
        % Integrate faces normal to the x direction.
        tangential_breaks_1 = VBreaks;
        tangential_breaks_2 = WBreaks;
        p1 = pv;
        p2 = pw;
        fixed_param = double(strcmp(faceName,'x+'));   % u = 0 or 1

        [gp1, gw1] = grule(10*p1 + 5);
        [gp2, gw2] = grule(10*p2 + 5);

        n_ele_1 = length(VBreaks) - 1;
        n_ele_2 = length(WBreaks) - 1;

        for e2 = 1:n_ele_2
            w_node = [tangential_breaks_2(e2), tangential_breaks_2(e2+1)];
            J2 = (w_node(2) - w_node(1))/2;

            for e1 = 1:n_ele_1
                v_node = [tangential_breaks_1(e1), tangential_breaks_1(e1+1)];
                J1 = (v_node(2) - v_node(1))/2;

                edge_dofs_plus = get_face_dofs_x(faceName, e1, e2, m, n, l, pv, pw, knotV, knotW, VBreaks, WBreaks);
                edge_dofs = [edge_dofs_plus, edge_dofs_minus];
                n_edge_dofs = length(edge_dofs);

                face_jump_Ae    = zeros(n_edge_dofs, n_edge_dofs);
                face_average_Ae = zeros(n_edge_dofs, n_edge_dofs);

                for i = 1:length(gp1)
                    v = ((v_node(2)-v_node(1))*gp1(i) + v_node(1) + v_node(2))/2;
                    for j = 1:length(gp2)
                        w = ((w_node(2)-w_node(1))*gp2(j) + w_node(1) + w_node(2))/2;

                        u = fixed_param;
                        [F, DF_plus] = NurbsVolume(ConPts_o, weights_o, ...
                            knotU_o, pu_o, u, knotV_o, pv_o, v, knotW_o, pw_o, w);

                        t1 = DF_plus(:,2);
                        t2 = DF_plus(:,3);
                        normal = cross(t1, t2);
                        if strcmp(faceName,'x-')
                            if dot(normal, [-1;0;0]) < 0
                                normal = -normal;
                            end
                        else
                            if dot(normal, [1;0;0]) < 0
                                normal = -normal;
                            end
                        end
                        ds = norm(normal);
                        normal = normal / ds;

                        Jacobi = J1 * J2 * gw1(i) * gw2(j) * ds;

                        [basis_plus, basis_grad_plus] = eval_face_basis_x( ...
                            faceName, u, v, w, knotU, knotV, knotW, pu, pv, pw, DF_plus);

                        for k = 1:n_pw_basis
                            basis_minus(k)        = exp(-1i*2*pi/L*pw_index(k,:)*F)/sqrt(Omega_vol);
                            basis_grad_minus(k,:) = (1i*2*pi/L) * pw_index(k,:) * basis_minus(k)';
                        end

                        edge_jump    = [basis_plus; -basis_minus];
                        face_jump_Ae = face_jump_Ae + edge_jump * edge_jump' * Jacobi;

                        edge_average    = [basis_grad_plus; basis_grad_minus] * normal / 2;
                        edge_average    = edge_average.';
                        face_average_Ae = face_average_Ae + edge_jump * edge_average * Jacobi;
                    end
                end

                P(edge_dofs, edge_dofs) = P(edge_dofs, edge_dofs) + face_jump_Ae;
                S(edge_dofs, edge_dofs) = S(edge_dofs, edge_dofs) + face_average_Ae;
            end
        end

    case {'y-','y+'}
        % Integrate faces normal to the y direction.
        tangential_breaks_1 = UBreaks;
        tangential_breaks_2 = WBreaks;
        p1 = pu;
        p2 = pw;
        fixed_param = double(strcmp(faceName,'y+'));   % v = 0 or 1

        [gp1, gw1] = grule(10*p1 + 5);
        [gp2, gw2] = grule(10*p2 + 5);

        n_ele_1 = length(UBreaks) - 1;
        n_ele_2 = length(WBreaks) - 1;

        for e2 = 1:n_ele_2
            w_node = [tangential_breaks_2(e2), tangential_breaks_2(e2+1)];
            J2 = (w_node(2) - w_node(1))/2;

            for e1 = 1:n_ele_1
                u_node = [tangential_breaks_1(e1), tangential_breaks_1(e1+1)];
                J1 = (u_node(2) - u_node(1))/2;

                edge_dofs_plus = get_face_dofs_y(faceName, e1, e2, m, n, l, pu, pw, knotU, knotW, UBreaks, WBreaks);
                edge_dofs = [edge_dofs_plus, edge_dofs_minus];
                n_edge_dofs = length(edge_dofs);

                face_jump_Ae    = zeros(n_edge_dofs, n_edge_dofs);
                face_average_Ae = zeros(n_edge_dofs, n_edge_dofs);

                for i = 1:length(gp1)
                    u = ((u_node(2)-u_node(1))*gp1(i) + u_node(1) + u_node(2))/2;
                    for j = 1:length(gp2)
                        w = ((w_node(2)-w_node(1))*gp2(j) + w_node(1) + w_node(2))/2;

                        v = fixed_param;
                        [F, DF_plus] = NurbsVolume(ConPts_o, weights_o, ...
                            knotU_o, pu_o, u, knotV_o, pv_o, v, knotW_o, pw_o, w);

                        t1 = DF_plus(:,1);
                        t2 = DF_plus(:,3);
                        normal = cross(t1, t2);
                        if strcmp(faceName,'y-')
                            if dot(normal, [0;-1;0]) < 0
                                normal = -normal;
                            end
                        else
                            if dot(normal, [0;1;0]) < 0
                                normal = -normal;
                            end
                        end
                        ds = norm(normal);
                        normal = normal / ds;

                        Jacobi = J1 * J2 * gw1(i) * gw2(j) * ds;

                        [basis_plus, basis_grad_plus] = eval_face_basis_y( ...
                            faceName, u, v, w, knotU, knotV, knotW, pu, pv, pw, DF_plus);

                        for k = 1:n_pw_basis
                            basis_minus(k)        = exp(-1i*2*pi/L*pw_index(k,:)*F)/sqrt(Omega_vol);
                            basis_grad_minus(k,:) = (1i*2*pi/L) * pw_index(k,:) * basis_minus(k)';
                        end

                        edge_jump    = [basis_plus; -basis_minus];
                        face_jump_Ae = face_jump_Ae + edge_jump * edge_jump' * Jacobi;

                        edge_average    = [basis_grad_plus; basis_grad_minus] * normal / 2;
                        edge_average    = edge_average.';
                        face_average_Ae = face_average_Ae + edge_jump * edge_average * Jacobi;
                    end
                end

                P(edge_dofs, edge_dofs) = P(edge_dofs, edge_dofs) + face_jump_Ae;
                S(edge_dofs, edge_dofs) = S(edge_dofs, edge_dofs) + face_average_Ae;
            end
        end

    case {'z-','z+'}
        % Integrate faces normal to the z direction.
        tangential_breaks_1 = UBreaks;
        tangential_breaks_2 = VBreaks;
        p1 = pu;
        p2 = pv;
        fixed_param = double(strcmp(faceName,'z+'));   % w = 0 or 1

        [gp1, gw1] = grule(10*p1 + 5);
        [gp2, gw2] = grule(10*p2 + 5);

        n_ele_1 = length(UBreaks) - 1;
        n_ele_2 = length(VBreaks) - 1;

        for e2 = 1:n_ele_2
            v_node = [tangential_breaks_2(e2), tangential_breaks_2(e2+1)];
            J2 = (v_node(2) - v_node(1))/2;

            for e1 = 1:n_ele_1
                u_node = [tangential_breaks_1(e1), tangential_breaks_1(e1+1)];
                J1 = (u_node(2) - u_node(1))/2;

                edge_dofs_plus = get_face_dofs_z(faceName, e1, e2, m, n, l, pu, pv, knotU, knotV, UBreaks, VBreaks);
                edge_dofs = [edge_dofs_plus, edge_dofs_minus];
                n_edge_dofs = length(edge_dofs);

                face_jump_Ae    = zeros(n_edge_dofs, n_edge_dofs);
                face_average_Ae = zeros(n_edge_dofs, n_edge_dofs);

                for i = 1:length(gp1)
                    u = ((u_node(2)-u_node(1))*gp1(i) + u_node(1) + u_node(2))/2;
                    for j = 1:length(gp2)
                        v = ((v_node(2)-v_node(1))*gp2(j) + v_node(1) + v_node(2))/2;

                        w = fixed_param;
                        [F, DF_plus] = NurbsVolume(ConPts_o, weights_o, ...
                            knotU_o, pu_o, u, knotV_o, pv_o, v, knotW_o, pw_o, w);

                        t1 = DF_plus(:,1);
                        t2 = DF_plus(:,2);
                        normal = cross(t1, t2);
                        if strcmp(faceName,'z-')
                            if dot(normal, [0;0;-1]) < 0
                                normal = -normal;
                            end
                        else
                            if dot(normal, [0;0;1]) < 0
                                normal = -normal;
                            end
                        end
                        ds = norm(normal);
                        normal = normal / ds;

                        Jacobi = J1 * J2 * gw1(i) * gw2(j) * ds;

                        [basis_plus, basis_grad_plus] = eval_face_basis_z( ...
                            faceName, u, v, w, knotU, knotV, knotW, pu, pv, pw, DF_plus);

                        for k = 1:n_pw_basis
                            basis_minus(k)        = exp(-1i*2*pi/L*pw_index(k,:)*F)/sqrt(Omega_vol);
                            basis_grad_minus(k,:) = (1i*2*pi/L) * pw_index(k,:) * basis_minus(k)';
                        end

                        edge_jump    = [basis_plus; -basis_minus];
                        face_jump_Ae = face_jump_Ae + edge_jump * edge_jump' * Jacobi;

                        edge_average    = [basis_grad_plus; basis_grad_minus] * normal / 2;
                        edge_average    = edge_average.';
                        face_average_Ae = face_average_Ae + edge_jump * edge_average * Jacobi;
                    end
                end

                P(edge_dofs, edge_dofs) = P(edge_dofs, edge_dofs) + face_jump_Ae;
                S(edge_dofs, edge_dofs) = S(edge_dofs, edge_dofs) + face_average_Ae;
            end
        end

    otherwise
        error('Unknown faceName = %s', faceName);
end

end

function edge_dofs_plus = get_face_dofs_x(faceName, e_v, e_w, m, n, ~, pv, pw, knotV, knotW, VBreaks, WBreaks)
% Collect control-point DOFs on an x-normal face.
jspan = findspan(knotV, pv, VBreaks(e_v));
kspan = findspan(knotW, pw, WBreaks(e_w));

j_list = (jspan-pv):jspan;
k_list = (kspan-pw):kspan;

if strcmp(faceName,'x-')
    i_list = [1, 2];
else
    i_list = [m-1, m];
end

edge_dofs_plus = zeros(1, 2*(pv+1)*(pw+1));
cnt = 0;
for kk = 1:(pw+1)
    for jj = 1:(pv+1)
        for ii = 1:2
            cnt = cnt + 1;
            edge_dofs_plus(cnt) = id3(i_list(ii), j_list(jj), k_list(kk), m, n);
        end
    end
end
end

function edge_dofs_plus = get_face_dofs_y(faceName, e_u, e_w, m, n, ~, pu, pw, knotU, knotW, UBreaks, WBreaks)
% Collect control-point DOFs on a y-normal face.
ispan = findspan(knotU, pu, UBreaks(e_u));
kspan = findspan(knotW, pw, WBreaks(e_w));

i_list = (ispan-pu):ispan;
k_list = (kspan-pw):kspan;

if strcmp(faceName,'y-')
    j_list = [1, 2];
else
    j_list = [n-1, n];
end

edge_dofs_plus = zeros(1, 2*(pu+1)*(pw+1));
cnt = 0;
for kk = 1:(pw+1)
    for ii = 1:(pu+1)
        for jj = 1:2
            cnt = cnt + 1;
            edge_dofs_plus(cnt) = id3(i_list(ii), j_list(jj), k_list(kk), m, n);
        end
    end
end
end

function edge_dofs_plus = get_face_dofs_z(faceName, e_u, e_v, m, n, l, pu, pv, knotU, knotV, UBreaks, VBreaks)
% Collect control-point DOFs on a z-normal face.
ispan = findspan(knotU, pu, UBreaks(e_u));
jspan = findspan(knotV, pv, VBreaks(e_v));

i_list = (ispan-pu):ispan;
j_list = (jspan-pv):jspan;

if strcmp(faceName,'z-')
    k_list = [1, 2];
else
    k_list = [l-1, l];
end

edge_dofs_plus = zeros(1, 2*(pu+1)*(pv+1));
cnt = 0;
for jj = 1:(pv+1)
    for ii = 1:(pu+1)
        for kk = 1:2
            cnt = cnt + 1;
            edge_dofs_plus(cnt) = id3(i_list(ii), j_list(jj), k_list(kk), m, n);
        end
    end
end
end

function [basis_plus, basis_grad_plus] = eval_face_basis_x(faceName, u, v, w, knotU, knotV, knotW, pu, pv, pw, DF_plus)
% Evaluate basis values and gradients on an x-normal face.
Uders = bspbasisDers(knotU, pu, u, 1);
Vders = bspbasisDers(knotV, pv, v, 1);
Wders = bspbasisDers(knotW, pw, w, 1);

if strcmp(faceName,'x-')
    Nu  = Uders(1,1:2)';
    DNu = Uders(2,1:2)';
else
    Nu  = Uders(1,end-1:end)';
    DNu = Uders(2,end-1:end)';
end

Nv  = Vders(1,:)';
DNv = Vders(2,:)';
Nw  = Wders(1,:)';
DNw = Wders(2,:)';

basis_plus = zeros(2*(pv+1)*(pw+1),1);
basis_grad_plus = zeros(2*(pv+1)*(pw+1),3);

cnt = 0;
for kk = 1:(pw+1)
    for jj = 1:(pv+1)
        for ii = 1:2
            cnt = cnt + 1;

            val_u = Nu(ii);
            val_v = Nv(jj);
            val_w = Nw(kk);

            d_u = DNu(ii) * val_v * val_w;
            d_v = val_u * DNv(jj) * val_w;
            d_w = val_u * val_v * DNw(kk);

            basis_plus(cnt) = val_u * val_v * val_w;
            basis_grad_plus(cnt,:) = [d_u, d_v, d_w] / DF_plus;
        end
    end
end
end

function [basis_plus, basis_grad_plus] = eval_face_basis_y(faceName, u, v, w, knotU, knotV, knotW, pu, pv, pw, DF_plus)
% Evaluate basis values and gradients on a y-normal face.
Uders = bspbasisDers(knotU, pu, u, 1);
Vders = bspbasisDers(knotV, pv, v, 1);
Wders = bspbasisDers(knotW, pw, w, 1);

Nu  = Uders(1,:)';
DNu = Uders(2,:)';

if strcmp(faceName,'y-')
    Nv  = Vders(1,1:2)';
    DNv = Vders(2,1:2)';
else
    Nv  = Vders(1,end-1:end)';
    DNv = Vders(2,end-1:end)';
end

Nw  = Wders(1,:)';
DNw = Wders(2,:)';

basis_plus = zeros(2*(pu+1)*(pw+1),1);
basis_grad_plus = zeros(2*(pu+1)*(pw+1),3);

cnt = 0;
for kk = 1:(pw+1)
    for ii = 1:(pu+1)
        for jj = 1:2
            cnt = cnt + 1;

            val_u = Nu(ii);
            val_v = Nv(jj);
            val_w = Nw(kk);

            d_u = DNu(ii) * val_v * val_w;
            d_v = val_u * DNv(jj) * val_w;
            d_w = val_u * val_v * DNw(kk);

            basis_plus(cnt) = val_u * val_v * val_w;
            basis_grad_plus(cnt,:) = [d_u, d_v, d_w] / DF_plus;
        end
    end
end
end

function [basis_plus, basis_grad_plus] = eval_face_basis_z(faceName, u, v, w, knotU, knotV, knotW, pu, pv, pw, DF_plus)
% Evaluate basis values and gradients on a z-normal face.
Uders = bspbasisDers(knotU, pu, u, 1);
Vders = bspbasisDers(knotV, pv, v, 1);
Wders = bspbasisDers(knotW, pw, w, 1);

Nu  = Uders(1,:)';
DNu = Uders(2,:)';

Nv  = Vders(1,:)';
DNv = Vders(2,:)';

if strcmp(faceName,'z-')
    Nw  = Wders(1,1:2)';
    DNw = Wders(2,1:2)';
else
    Nw  = Wders(1,end-1:end)';
    DNw = Wders(2,end-1:end)';
end

basis_plus = zeros(2*(pu+1)*(pv+1),1);
basis_grad_plus = zeros(2*(pu+1)*(pv+1),3);

cnt = 0;
for jj = 1:(pv+1)
    for ii = 1:(pu+1)
        for kk = 1:2
            cnt = cnt + 1;

            val_u = Nu(ii);
            val_v = Nv(jj);
            val_w = Nw(kk);

            d_u = DNu(ii) * val_v * val_w;
            d_v = val_u * DNv(jj) * val_w;
            d_w = val_u * val_v * DNw(kk);

            basis_plus(cnt) = val_u * val_v * val_w;
            basis_grad_plus(cnt,:) = [d_u, d_v, d_w] / DF_plus;
        end
    end
end
end

function idx = id3(i,j,k,m,n)
% Map three tensor indices to one linear index.
idx = i + (j-1)*m + (k-1)*m*n;
end
