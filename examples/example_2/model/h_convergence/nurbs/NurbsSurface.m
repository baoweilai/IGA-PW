function [S,DF,W,DWu,DWv]=NurbsSurface(ConPts,wights,knotU,pu,u,knotV,pv,v)
% Evaluate a NURBS surface.

[~,~,ndim]=size(ConPts);
Pw=WightedConPtsSurface(ConPts,wights);
Sw=PointOnBspSurface(Pw,knotU,pu,u,knotV,pv,v);
S=Project(Sw);
W=Sw(end);
[DSwu,DSwv]=bspSurfaceder(Pw,knotU,pu,u,knotV,pv,v);
DAu=DSwu(1:ndim);DAv=DSwv(1:ndim);DWu=DSwu(end);DWv=DSwv(end);
DSu=(DAu-DWu*S)/W;
DSv=(DAv-DWv*S)/W;
DF=[DSu,DSv];
