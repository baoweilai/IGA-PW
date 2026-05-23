function S=PointOnNurbsSurface(ConPts,wights,knotU,pu,u,knotV,pv,v)
%Evaluate a point on a NURBS surface.


[~,~,ndim]=size(ConPts);
Pw=WightedConPtsSurface(ConPts,wights);
Sw=PointOnBspSurface(Pw,knotU,pu,u,knotV,pv,v);
S=Project(Sw);
