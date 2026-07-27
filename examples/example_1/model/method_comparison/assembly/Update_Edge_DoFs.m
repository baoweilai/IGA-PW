function nurbs_refine = Update_Edge_DoFs(n_dofs_1,nurbs_refine)
% Update edge degrees of freedom.

% From Omega_1 to Omega_2, we have to update the index of dofs on Omega_2

nurbs_refine.bottom_edge_dofs =  n_dofs_1 +  nurbs_refine.bottom_edge_dofs  ;
nurbs_refine.right_edge_dofs  =  n_dofs_1 +  nurbs_refine.right_edge_dofs   ;
nurbs_refine.top_edge_dofs    =  n_dofs_1 +  nurbs_refine.top_edge_dofs     ;
nurbs_refine.left_edge_dofs   =  n_dofs_1 +  nurbs_refine.left_edge_dofs    ;

nurbs_refine.bottom_dofs = n_dofs_1 +  nurbs_refine.bottom_dofs;
nurbs_refine.top_dofs    = n_dofs_1 +  nurbs_refine.top_dofs;
nurbs_refine.left_dofs   = n_dofs_1 +  nurbs_refine.left_dofs;
nurbs_refine.right_dofs  = n_dofs_1 +  nurbs_refine.right_dofs;

end
