function L2_err   = Compute_L2_Error(nurbs_original,nurbs_refine,uh,ref_nurbs,ref_uh)
% Compute the L2 error on the comparison grid.

% This script is used to compute the $L^2$ norm, $H^1$ semi-norm and $H^2$
% semi-norm  error between $u$ and the IGA solution $u_h$.

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU;
knotV_o    = nurbs_original.knotV;
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;

pu   =  nurbs_refine.pu;
pv   =  nurbs_refine.pv;
Element=nurbs_refine.Element;
Coordinate = nurbs_refine.Coordinate;
Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
NoEs=nurbs_refine.NoEs;

u_np = pu + 5 ;   % It seems that if pu >=3, u_np = pu is enough!!!
v_np = pv + 5;

[gp_u,gw_u]=grule(u_np); % The quadrature points and quadrature weights in the [-1,1].

gp_v = gp_u;
gw_v = gw_u;

L2_err = 0;

for e=1:NoEs % Loop for each element;
    row=Element(e,:);
    ue=Coordinate(e,1:2);
    ve=Coordinate(e,3:4);
    U_ij = uh(row);
    U_ij = U_ij';
    for j=1:v_np
        vJ = (ve(2) - ve(1))/2;
        v = vJ*gp_v(j) + (ve(1)+ve(2))/2;
        Nv = bsplinebasis(Vbar,pv,v); % A column vector
        for i=1:u_np
            uJ = (ue(2) - ue(1))/2;
            u = uJ*gp_u(i)+(ue(1)+ue(2))/2;
            Nu = bsplinebasis(Ubar,pu,u); % A column vector
            [~,~, ~,~,DF]=NurbsSurfaceDers(ConPts_o,knotU_o,knotV_o,weights_o,pu_o,u,pv_o,v);
            J=abs(det(DF))*gw_u(i)*uJ*gw_v(j)*vJ;
            B= Nu*(Nv');
            B= B(:);

            ref_uh_at_uv   = reference_uh(ref_nurbs,ref_uh,u,v);
            L2_err = L2_err + abs( ref_uh_at_uv - U_ij*B ).^2*J;
        end
    end
end
L2_err  = sqrt(L2_err );
end
