function [P,S] =  IGA_DG_Top_Bottom_Edge_Assemble(nurbs_original_1,nurbs_original_2,nurbs_refine_1,nurbs_refine_2,n_dofs)
%Assemble matrices or interface terms for the method.

% Now the two sub-domains are the upper one (v=0)  and lower one (v=1), the
% interface is at v = 0 for upper one, and v = 1 for the upper one


ConPts_o_1   =  nurbs_original_1.ConPts;
weights_o_1  =  nurbs_original_1.weights;
knotU_o_1    =  nurbs_original_1.knotU;
knotV_o_1    =  nurbs_original_1.knotV;
pu_o_1       =  nurbs_original_1.pu;
pv_o_1       =  nurbs_original_1.pv;

knotU_1 = nurbs_refine_1.Ubar;
knotV_1 = nurbs_refine_1.Vbar;
pu_1    = nurbs_refine_1.pu;
pv_1    = nurbs_refine_1.pv;

ConPts_o_2   =  nurbs_original_2.ConPts;
weights_o_2  =  nurbs_original_2.weights;
knotU_o_2    =  nurbs_original_2.knotU;
knotV_o_2    =  nurbs_original_2.knotV;
pu_o_2       =  nurbs_original_2.pu;
pv_o_2       =  nurbs_original_2.pv;

knotU_2 = nurbs_refine_2.Ubar;
knotV_2 = nurbs_refine_2.Vbar;
pu_2    = nurbs_refine_2.pu;
pv_2    = nurbs_refine_2.pv;

S = sparse(n_dofs,n_dofs);
P = S;
% The plus sign points from the lower domain to upper element
UBreaks_1 = nurbs_refine_1.UBreaks;
UBreaks_2 = nurbs_refine_2.UBreaks;

top_edge_dofs_1 = nurbs_refine_1.top_edge_dofs;
bottom_edge_dofs_2  = nurbs_refine_2.bottom_edge_dofs;

UBreaks_merge = [UBreaks_1, UBreaks_2];
UBreaks_merge = unique(UBreaks_merge);
uNoEs = length(UBreaks_merge) - 1;

bottom_edge_node_1 = zeros(uNoEs,2);
for i=1:uNoEs
    bottom_edge_node_1(i,:) = [UBreaks_merge(i),UBreaks_merge(i+1)];
end


edge_element_1 =  Merge_edge_elements_idx(UBreaks_1,UBreaks_merge);
edge_element_2 =  Merge_edge_elements_idx(UBreaks_2,UBreaks_merge);

[gp,gw] = grule(pu_1+3); % The Gaussian quadrature rule on [-1,1].

n_gp = length(gp);

% As the interface is v = 1, so we consider v = 1
v_top_1 = 1;
v_bottom_2 = 0;

for e=1:uNoEs

    edge_idx_plus  = edge_element_1(e);
    edge_idx_minus = edge_element_2(e);

    edge_dofs_plus  = top_edge_dofs_1(edge_idx_plus,:);
    edge_dofs_minus = bottom_edge_dofs_2(edge_idx_minus,:);

    edge_dofs = [edge_dofs_plus,edge_dofs_minus];

    n_edge_dofs = length(edge_dofs);

    ue = bottom_edge_node_1(e,:);
    a  = ue(1);    b = ue(2);
    J1 = (b-a)/2;  % The Jacobian from [-1,1] to [a,b].
    edge_jump_Ae    = zeros(n_edge_dofs,n_edge_dofs);
    edge_average_Ae = zeros(n_edge_dofs,n_edge_dofs);

    for i=1:n_gp
        u  = ((b-a)*gp(i) +a+b)/2;
        [F,DF_plus]   =  NurbsSurface(ConPts_o_1 ,weights_o_1,knotU_o_1 ,pu_o_1,u,knotV_o_1,pv_o_1,v_top_1);% Left patch

        tau    = DF_plus(:,1);
        ds     = norm(tau);
        normal = [-tau(2);tau(1)]/ds;


        Jacobi = J1*gw(i)*ds; % The interface is u = 1, now along the v-direction

        Uders_plus  = bspbasisDers(knotU_1,pu_1,u,1);
        Nu_plus = Uders_plus(1,:)';     DNu_plus = Uders_plus(2,:)';
        Vders_plus  = bspbasisDers(knotV_1,pv_1,v_top_1,1);
        Nv_plus = Vders_plus(1,end-1:end);       DNv_plus = Vders_plus(2,end-1:end);

        basis_plus = Nu_plus*Nv_plus;    basis_plus = basis_plus(:);
        DNu_v_plus = DNu_plus*Nv_plus;   DNu_v_plus = DNu_v_plus(:);
        DNv_u_plus = Nu_plus*DNv_plus;   DNv_u_plus = DNv_u_plus(:);
        basis_grad_plus = [DNu_v_plus,DNv_u_plus]/DF_plus;


        Uders_minus  = bspbasisDers(knotU_2,pu_2,u,1); % For Omega 2, it is the left boundary, u = 0.
        Nu_minus     = Uders_minus(1,:)';     DNu_minus = Uders_minus(2,:)';
        Vders_minus  = bspbasisDers(knotV_2,pv_2,v_bottom_2,1);
        Nv_minus     = Vders_minus(1,1:2);              DNv_minus = Vders_minus(2,1:2);

        [~,DF_minus] =  NurbsSurface(ConPts_o_2 ,weights_o_2,knotU_o_2 ,pu_o_2, u ,knotV_o_2,pv_o_2 ,v_bottom_2);% Right patch
        basis_minus = Nu_minus*Nv_minus;    basis_minus = basis_minus(:);
        DNu_v_minus = DNu_minus*Nv_minus;   DNu_v_minus = DNu_v_minus(:);
        DNv_u_minus = Nu_minus*DNv_minus;   DNv_u_minus = DNv_u_minus(:);
        basis_grad_minus = [DNu_v_minus,DNv_u_minus]/DF_minus;

        edge_jump    = [basis_plus; - basis_minus];
        edge_jump_Ae = edge_jump_Ae + edge_jump*edge_jump'*Jacobi;


        edge_average    = [basis_grad_plus; basis_grad_minus]*normal/2; % Now it is a column vector
        edge_average    = edge_average'; % Now it is a row vector
        edge_average_Ae = edge_average_Ae + edge_jump*edge_average*Jacobi;


    end


    P(edge_dofs,edge_dofs) =  P(edge_dofs,edge_dofs) + edge_jump_Ae;

    S(edge_dofs,edge_dofs) =  S(edge_dofs,edge_dofs) + edge_average_Ae;

end


end
