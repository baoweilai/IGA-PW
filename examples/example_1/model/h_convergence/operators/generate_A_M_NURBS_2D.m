function [A,M] = generate_A_M_NURBS_2D(nurbs_original,nurbs_refine,p_Vr,n_pw_Vr, L,n_gp, Example)
%Assemble 2-D NURBS matrices.

ConPts_o     = nurbs_original.ConPts ;
weights_o    = nurbs_original.weights;
knotU_o      = nurbs_original.knotU;
knotV_o      = nurbs_original.knotV;
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;
Element=nurbs_refine.Element;
knotU=nurbs_refine.Ubar;
knotV=nurbs_refine.Vbar;
Coordinate = nurbs_refine.Coordinate;
Fhat = @(x,a,b) ( (b-a)*x+a+b )/2;
NoEs=nurbs_refine.NoEs;
n_dofs=nurbs_refine.n_dofs;

pu = nurbs_refine.pu;
pv = nurbs_refine.pv;
[gp,gw] = grule(n_gp);
n_ele_dofs = (pu+1)*(pv+1);

Ae =zeros(n_ele_dofs,n_ele_dofs); % Element stiffness matrix
Me =zeros(n_ele_dofs,n_ele_dofs); % Element mass matrix

A_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
M_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
row_index = A_value;
column_index = A_value;
global_index = 1;

for e = 1:NoEs
    Ae = 0*Ae;  Me = 0*Me; % Fe =0*Fe;     Me_ext = 0*Me_ext;
    ue = Coordinate(e,1:2);   J1 = (ue(2) - ue(1) )/2;
    ve = Coordinate(e,3:4);   J2 = (ve(2) - ve(1) )/2;
    row = Element(e,:);
    for i=1:n_gp
        u       = Fhat(gp(i),ue(1),ue(2));
        uders   = bspbasisDers(knotU,pu,u,1);
        Nu  = uders(1,:)';
        DNu = uders(2,:)';
        for j=1:n_gp
            v       = Fhat(gp(j),ve(1),ve(2));
            vders   = bspbasisDers(knotV,pv,v,1);
            Nv  = vders(1,:);
            DNv = vders(2,:);
            [F,DF] = NurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v);
            basis_funcs = Nu*Nv;  basis_funcs = basis_funcs(:);
            DNu_v = DNu*Nv; DNu_v = DNu_v(:);
            DNv_u = Nu*DNv; DNv_u = DNv_u(:);
            basis_grad = [DNu_v,DNv_u]/DF;
            Jacobian = J1*J2*abs(det(DF))*gw(i)*gw(j);
            Mass_mat_ele = basis_funcs*basis_funcs'*Jacobian;
            if strcmp( Example, 'Example_1')
                V_ext = Vr_2D_Example_1(p_Vr,L,n_pw_Vr,F(1),F(2));
            else
                V_ext = Vr_2D_Example_2(p_Vr,L,n_pw_Vr,F(1),F(2));
            end
            Ae = Ae + basis_grad*basis_grad'*Jacobian/2 + V_ext*Mass_mat_ele;
            Me = Me + Mass_mat_ele;
        end
    end

    for i1=1:n_ele_dofs
        for j1=1:n_ele_dofs
            row_index(global_index) = row(i1);
            column_index(global_index) = row(j1);
            A_value(global_index) =   Ae(i1,j1);
            M_value(global_index) =   Me(i1,j1);
            global_index = global_index + 1;
        end
    end
end
A =sparse(row_index,column_index,A_value,n_dofs,n_dofs);
M =sparse(row_index,column_index,M_value,n_dofs,n_dofs);
end
