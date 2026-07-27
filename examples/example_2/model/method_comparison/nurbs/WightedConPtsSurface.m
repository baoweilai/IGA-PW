function Pw=WightedConPtsSurface(P,w)
% Convert surface control points to homogeneous form.

[m,n,ndim]=size(P);
Pw=zeros(m,n,ndim+1);
Pw(:,:,end)=w;
for i=1:ndim
    Pw(:,:,i)=P(:,:,i).*w;
end
end
