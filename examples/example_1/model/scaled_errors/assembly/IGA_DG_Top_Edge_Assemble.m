function [P,S] =  IGA_DG_Top_Edge_Assemble(nurbs_original, nurbs_refine, pw_index,plane_wave_dofs_index, L,  n_dofs)
%Assemble matrices or interface terms for the method.

% Now the two sub-domains are the upper one (v=0)  and lower one (v=1), the
% interface is at v = 0 for upper one, and v = 1 for the upper one

DIM = 2;


ConPts_o   =  nurbs_original.ConPts;
weights_o  =  nurbs_original.weights;
knotU_o    =  nurbs_original.knotU;
knotV_o    =  nurbs_original.knotV;
pu_o       =  nurbs_original.pu;
pv_o       =  nurbs_original.pv;


knotU     = nurbs_refine.Ubar;
knotV     = nurbs_refine.Vbar;
pu        = nurbs_refine.pu;
pv        = nurbs_refine.pv;


S = sparse(n_dofs,n_dofs);

P = S;

% The plus sign points from the lower domain (inner domain) to upper (outer) domain

UBreaks       = nurbs_refine.UBreaks;
uNoEs         = length(UBreaks) - 1;

top_edge_dofs = nurbs_refine.top_edge_dofs;


top_edge_node = zeros(uNoEs,2);

for i=1:uNoEs
    top_edge_node(i,:) = [UBreaks(i),UBreaks(i+1)];
end


[gp,gw] = grule(10*pu+5); % The Gaussian quadrature rule on [-1,1].

n_gp = length(gp);


% As the interface is top boundary of inner domain, then v = 1

v_top_inner = 1;

%% For plane wave functions
edge_dofs_minus  =  plane_wave_dofs_index; % To be done
n_pw_basis       =  size(pw_index,1);
basis_grad_minus =  zeros(n_pw_basis,DIM);
basis_minus      =  zeros(n_pw_basis,1);
Omega_area       =  L*L;


for e=1:uNoEs

    edge_dofs_plus  = top_edge_dofs(e,:);
    edge_dofs       = [edge_dofs_plus,edge_dofs_minus];

    n_edge_dofs = length(edge_dofs);

    ue = top_edge_node(e,:);
    a  = ue(1);    b = ue(2);
    J1 = (b-a)/2;  % The Jacobian from [-1,1] to [a,b].
    edge_jump_Ae    = zeros(n_edge_dofs,n_edge_dofs);
    edge_average_Ae = zeros(n_edge_dofs,n_edge_dofs);

    for i=1:n_gp
        u  = ((b-a)*gp(i) +a+b)/2;
        [F,DF_plus]   =  NurbsSurface(ConPts_o ,weights_o,knotU_o ,pu_o,u,knotV_o,pv_o,v_top_inner);% Left patch

        tau    = DF_plus(:,1);
        ds     = norm(tau);
        normal = [-tau(2);tau(1)]/ds;


        Jacobi = J1*gw(i)*ds; % The interface is u = 1, now along the v-direction

        Uders_plus  = bspbasisDers(knotU,pu,u,1);
        Nu_plus = Uders_plus(1,:)';              DNu_plus = Uders_plus(2,:)';
        Vders_plus  = bspbasisDers(knotV,pv,v_top_inner,1);
        Nv_plus = Vders_plus(1,end-1:end);       DNv_plus = Vders_plus(2,end-1:end);

        basis_plus = Nu_plus*Nv_plus;    basis_plus = basis_plus(:);
        DNu_v_plus = DNu_plus*Nv_plus;   DNu_v_plus = DNu_v_plus(:);
        DNv_u_plus = Nu_plus*DNv_plus;   DNv_u_plus = DNv_u_plus(:);
        basis_grad_plus = [DNu_v_plus,DNv_u_plus]/DF_plus;


        %% The basis functions in outer domain are plane wave functions, which do not vanish on the boundary

        for k = 1:n_pw_basis

            basis_minus(k)        = exp(-1i*2*pi/L*pw_index(k,:)*F)/sqrt(Omega_area);
            basis_grad_minus(k,:) = (1i*2*pi/L)*pw_index(k,:)*basis_minus(k)';

        end


        edge_jump    = [basis_plus; - basis_minus];
        edge_jump_Ae = edge_jump_Ae + edge_jump*edge_jump'*Jacobi;


        edge_average    = [basis_grad_plus; basis_grad_minus]*normal/2; % Now it is a column vector
        edge_average    = edge_average.';
        edge_average_Ae = edge_average_Ae + edge_jump*edge_average*Jacobi;


    end


    P(edge_dofs,edge_dofs) =  P(edge_dofs,edge_dofs) + edge_jump_Ae;

    S(edge_dofs,edge_dofs) =  S(edge_dofs,edge_dofs) + edge_average_Ae;

end


end
