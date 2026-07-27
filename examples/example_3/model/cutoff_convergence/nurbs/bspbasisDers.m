function ders=bspbasisDers(U,p,u,n)
% Evaluate B-spline basis derivatives.

ders=zeros(n+1,p+1);
ndu=AllBasisFuns(U,p,u);
ders(1,:)= ndu(:,p+1)';
a=zeros(2,p+1);

for r=0:p
    s1=1; s2=2; a(1,1)=1.0;
    for k=1:n
        d=0;
        rk=r-k;
        pk=p-k;

        if(rk>=0)
            a(s2,1)=a(s1,1)/ndu(pk+2,rk+1);
            d=d+a(s2,1)*ndu(rk+1,pk+1);
            j1=1;
        else
            j1=-rk;
        end

        if(r<=pk)
            a(s2,k+1)=-a(s1,k)/ndu(pk+2,r+1);
            d=d+a(s2,k+1)*ndu(r+1,pk+1);j2=k-1;
        else
            j2=p-r;
        end

        for j=j1:j2
            a(s2,j+1)=(a(s1,j+1)-a(s1,j))/ndu(pk+2,rk+j+1);
            d=d+a(s2,j+1)*ndu(rk+j+1,pk+1);
        end

        ders(k+1,r+1)=d;
        j=s1;s1=s2;s2=j;
    end
end

r=p;
for k=1:n
    ders(k+1,:)=r*ders(k+1,:);   r=r*(p-k);
end
