function S=PointOnBspSurface(P,U,p,u,V,q,v)
%Evaluate a point on a B-spline surface.
ndim=size(P,3);
uspan=findspan(U,p,u);Nu=bsplinebasis(U,p,u);
vspan=findspan(V,q,v);Nv=bsplinebasis(V,q,v);
S=zeros(ndim,1);

for i=1:ndim
  temp=reshape(P(uspan-p:uspan,vspan-q:vspan,i),p+1,q+1);
 S(i)=Nu'*temp*Nv;
end


%%==========================Test==================
% P(:,:,1)=[0 3 6 9;0 3 6 9;0 3 6 9]';
% P(:,:,2)=[0 0 0 0;2 2 2 2;4 4 4 4]';
% P(:,:,3)=[0 3 3 0;2 5 5 2;0 3 3 0]';

%% Section
%% Section

 % The Nurbs Book P116 3.8.

%% Section

%% S(i)=Nu'*temp*Nv;
