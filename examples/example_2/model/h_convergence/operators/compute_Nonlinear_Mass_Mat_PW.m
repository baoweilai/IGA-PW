function H = compute_Nonlinear_Mass_Mat_PW(L, Nc, inner_domains,rho, u0_h)
%Assemble the plane-wave nonlinear mass matrix.


%============================= initiate the data =========================%
DIM = 2;
n = 0; % number of basis satisfying |k|<=Nc
N = floor(Nc);

p = zeros((2*N+1)^2,DIM);

for ii = -N:N
    m = sqrt(N^2 - ii^2);
    m = floor(m);
    for j=-m:m
        n=n+1;
        p(n,:)=[ii,j];

    end
end

p = p(1:n,:);


% evaluate the discrete potential


dx = 1e-1;
m  = L/dx;
m  = floor(m);
dy = dx;
vext = zeros(m,m);


for ii = 1:m
    for j = 1:m

        x = -L/2 +  dx/2  + (ii-1)*L/m ;
        y = -L/2 +  dy/2  + (j-1)*L/m  ;
        vext(ii,j) = V_hydrogen(x,y,inner_domains,rho,u0_h,L,p, n);

    end
end


% inverse FFT for vext

vext_fft = ifftn(vext);


%======================== calculate the matrix element ===================%

H = zeros(n,n);

for ii=1:n
    for j=1:n
        dk  =  p(j,:) - p(ii,:);
        q_p = dk;


        % external potential H_{pq}=vext_fft(p-q)

        for k = 1:DIM
            if dk(k) < 0
                dk(k) = dk(k) + m;
            end

            dk(k) = dk(k) + 1;    % this is an index of fft matrix

            if dk(k)>m
                dk(k) = dk(k) - m;
            end

        end

        %          H(ii,j) = H(ii,j) + vext_fft(dk(1),dk(2)) *cos( q_p(1)*pi )*cos( q_p(2)*pi );    %% By using V(r) = 1/r;

        H(ii,j) = H(ii,j) + vext_fft( dk(1),dk(2) )*exp(1i* pi* q_p(1)* (1/m + 1) )* exp(1i* pi* q_p(2)* (1/m + 1) );

        %            H(ii,j) = H(ii,j) + M(ii,j);  % For V(1) = 1
    end
end

end


function Vr = V_hydrogen(x,y,inner_domains,rho,u0h,L,pw_index, n_pw_basis )
%Evaluate the hydrogen potential.

%                  a2 b2 c2 d2;   % The second inner sub-domain
%                  a3 b3 c3 d3;   % The third inner sub-domain
%                  an bn cn dn]   % The n-th inner sub-domain

n_inner_domains = size(inner_domains,1);

Omega_in_flag = false;

for i = 1:n_inner_domains

    x_a = inner_domains(i,1);
    x_b = inner_domains(i,2);
    y_c = inner_domains(i,3);
    y_d = inner_domains(i,4);

    if (x>x_a && x<x_b) && ( y>y_c && y<y_d) % The inner domain is [a,b]*[c,d]
        Omega_in_flag = true;
        break;
    end

end


if Omega_in_flag == true
    Vr = 0;
else

    uh_xy = 0;

    for k = 1:n_pw_basis

        uh_xy  = uh_xy + u0h(k)* exp( 1i * 2*pi*pw_index(k,:)*[x;y]/L )/sqrt(L*L);

    end

    Vr = rho( uh_xy );

end

end
