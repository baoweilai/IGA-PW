function [P,S] =  IGA_DG_Left_Right_Edge_Assemble(nurbs_original_1,nurbs_original_2,nurbs_refine_1,nurbs_refine_2,n_dofs)
% Assemble jump and average terms on a vertical patch interface.


% Read the geometry and refined-basis data for both patches.
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


% Merge both edge partitions and identify neighboring elements.
right_edge_dofs_1 = nurbs_refine_1.right_edge_dofs;
left_edge_dofs_2  = nurbs_refine_2.left_edge_dofs;


VBreaks_1 = nurbs_refine_1.VBreaks;
VBreaks_2 = nurbs_refine_2.VBreaks;


VBreaks_merge = [VBreaks_1, VBreaks_2];
VBreaks_merge = unique(VBreaks_merge);
vNoEs = length(VBreaks_merge) - 1;


right_edge_node_1 = zeros(vNoEs,2);

for i=1:vNoEs
right_edge_node_1(i,:) = [VBreaks_merge(i),VBreaks_merge(i+1)];
end


edge_element_1 =  Merge_edge_elements_idx(VBreaks_1,VBreaks_merge);
edge_element_2 =  Merge_edge_elements_idx(VBreaks_2,VBreaks_merge);


[gp,gw] = grule(pu_1+3); % The Gaussian quadrature rule on [-1,1].

n_gp = length(gp);


% Integrate the interface at u=1 on the first patch.
for e=1:vNoEs

    edge_idx_plus  = edge_element_1(e);
    edge_idx_minus = edge_element_2(e);

    edge_dofs_plus  = right_edge_dofs_1(edge_idx_plus,:);
    edge_dofs_minus = left_edge_dofs_2(edge_idx_minus,:);

    edge_dofs = [edge_dofs_plus,edge_dofs_minus];

    n_edge_dofs = length(edge_dofs);

    ve = right_edge_node_1(e,:);
    a  = ve(1);    b = ve(2);
    J1 = (b-a)/2;  % The Jacobian from [-1,1] to [a,b].
    edge_jump_Ae    = zeros(n_edge_dofs,n_edge_dofs);
    edge_average_Ae = zeros(n_edge_dofs,n_edge_dofs);

    for i=1:n_gp
    v  = ((b-a)*gp(i) +a+b)/2;
    [~,DF_plus]   =  NurbsSurface(ConPts_o_1 ,weights_o_1,knotU_o_1 ,pu_o_1,1,knotV_o_1,pv_o_1 ,v);% Left patch

    tau    = DF_plus(:,2);
    ds     = norm(tau);
    normal = [tau(2);-tau(1)]/ds;


    Jacobi = J1*gw(i)*ds; % Scale quadrature along the v direction.

    % Evaluate traces and physical gradients on both sides.
    Uders_plus  = bspbasisDers(knotU_1,pu_1,1,1);
    Nu_plus = Uders_plus(1,end-1:end)';     DNu_plus = Uders_plus(2,end-1:end)';
    Vders_plus  = bspbasisDers(knotV_1,pv_1,v,1);
    Nv_plus = Vders_plus(1,:);       DNv_plus = Vders_plus(2,:);

    basis_plus = Nu_plus*Nv_plus;    basis_plus = basis_plus(:);
    DNu_v_plus = DNu_plus*Nv_plus;   DNu_v_plus = DNu_v_plus(:);
    DNv_u_plus = Nu_plus*DNv_plus;   DNv_u_plus = DNv_u_plus(:);
    basis_grad_plus = [DNu_v_plus,DNv_u_plus]/DF_plus;


    Uders_minus  = bspbasisDers(knotU_2,pu_2,0,1); % The second patch trace uses u=0.
    Nu_minus     = Uders_minus(1,1:2)';     DNu_minus = Uders_minus(2,1:2)';
    Vders_minus  = bspbasisDers(knotV_2,pv_2,v,1);
    Nv_minus     = Vders_minus(1,:);              DNv_minus = Vders_minus(2,:);

    [~,DF_minus] =  NurbsSurface(ConPts_o_2 ,weights_o_2,knotU_o_2 ,pu_o_2,0,knotV_o_2,pv_o_2 ,v);% Right patch
    basis_minus = Nu_minus*Nv_minus;    basis_minus = basis_minus(:);
    DNu_v_minus = DNu_minus*Nv_minus;   DNu_v_minus = DNu_v_minus(:);
    DNv_u_minus = Nu_minus*DNv_minus;   DNv_u_minus = DNv_u_minus(:);
    basis_grad_minus = [DNu_v_minus,DNv_u_minus]/DF_minus;

    edge_jump    = [basis_plus; - basis_minus];


    edge_jump_Ae = edge_jump_Ae + edge_jump*edge_jump'*Jacobi;


    edge_average    = [basis_grad_plus; basis_grad_minus]*normal/2;
    edge_average    = edge_average'; % Form the averaged normal derivative as a row vector.
    edge_average_Ae = edge_average_Ae + edge_jump*edge_average*Jacobi;


    end


    % Add the local interface blocks to the global matrices.
    P(edge_dofs,edge_dofs) =  P(edge_dofs,edge_dofs) + edge_jump_Ae;

    S(edge_dofs,edge_dofs) =  S(edge_dofs,edge_dofs) + edge_average_Ae;

end


end
