function S=PointOnBspSurface(P,U,p,u,V,q,v)
% Evaluate a point on a B-spline surface.
ndim=size(P,3);
uspan=findspan(U,p,u);Nu=bsplinebasis(U,p,u);
vspan=findspan(V,q,v);Nv=bsplinebasis(V,q,v);
S=zeros(ndim,1);

% Contract the active control net with the tensor-product basis.
for i=1:ndim
  temp=reshape(P(uspan-p:uspan,vspan-q:vspan,i),p+1,q+1);
 S(i)=Nu'*temp*Nv;
end
