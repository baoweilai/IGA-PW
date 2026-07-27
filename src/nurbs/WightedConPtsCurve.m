function Pw=WightedConPtsCurve(P,w)
% Convert curve control points to homogeneous form.

[ndim,n]=size(P);
Pw=zeros(ndim+1,n);
Pw(end,:)=w;
for i=1:ndim
    Pw(i,:)=P(i,:).*w;
end

end
