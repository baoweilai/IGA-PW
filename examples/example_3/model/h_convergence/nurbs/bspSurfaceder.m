function [DSu,DSv]=bspSurfaceder(P,U,p,u,V,q,v)
% Evaluate B-spline surface derivatives.

uspan=findspan(U,p,u);
vspan=findspan(V,q,v);
ndim=size(P,3);
temp=P(uspan-p:uspan,vspan-q:vspan,:);
ConPtsUbar=zeros(p,q+1,ndim);
ConPtsVbar=zeros(p+1,q,ndim);
tempU=(U(uspan+1:uspan+p)-U(uspan-p+1:uspan))';
tempU=tempU*ones(1,q+1);
tempV=V(vspan+1:vspan+q)-V(vspan-q+1:vspan);
tempV=ones(p+1,1)*tempV;
for  i=1:ndim
    ConPtsUbar(:,:,i)=p*(temp(2:end,:,i)-temp(1:end-1,:,i))./tempU;
    ConPtsVbar(:,:,i)=q*(temp(:,2:end,i)-temp(:,1:end-1,i))./tempV;
end

% Evaluate the derivative in the u direction.
Ubar=U;
Ubar([1,end])=[];
Nu=bsplinebasis(Ubar,p-1,u);
Nv=bsplinebasis(V,q,v);
DSu=zeros(ndim,1);
for i=1:ndim
    temp=reshape(ConPtsUbar(:,:,i),p,q+1);
    DSu(i)=Nu'*temp*Nv;
end

% Evaluate the derivative in the v direction.
Vbar=V;
Vbar([1,end])=[];
Nu=bsplinebasis(U,p,u);
Nv=bsplinebasis(Vbar,q-1,v);
DSv=zeros(ndim,1);
for i=1:ndim
    temp=reshape(ConPtsVbar(:,:,i),p+1,q);
    DSv(i)=Nu'*temp*Nv;
end

end
