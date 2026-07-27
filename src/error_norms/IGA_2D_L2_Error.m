function L2_err =  IGA_2D_L2_Error( Refinement, Example, t, Nc, n_eigenvalues  )
% Compute the 2-D IGA L2 error.

% Build the circular plane-wave index set.
DIM = 2;
n_pw_basis = 0; % number of basis satisfying |k|<=Nc
N=floor(Nc);
p = zeros((2*N+1)^2,DIM);
for ii = -N:N
    m = sqrt(N^2 - ii^2);
    m = floor(m);
    for j=-m:m
        n_pw_basis=n_pw_basis+1;
        p(n_pw_basis,:)=[ii,j];
    end
end
p = p(1:n_pw_basis,:);

format long
DIM  = 2;

%% For the domain with one singularity
if strcmp(Example, 'Example_1') || strcmp(Example, 'Example_3')
    %% For the domain with one singularity
    L  = 4;    % The domain of Example 1 is Omega = [-L/2,L/2]^2
    a  = 0.2;  % The inner domain is [-a,a]^2
    n_nurbs_patch = 1;
    x  = [-a,a];
    y  = [-a,a];
    inner_domains_coordinates = [x,y];
end

%%  For the domain with two singularities: Omega = [-1,1]^2
if strcmp(Example, 'Example_2')
    L  = 4;       % The domain is [-L/2,L/2]
    x1 = -1;     x2 = 1;
    y1 =  0;     y2 = 0;
    n_nurbs_patch = 2;
    a = 0.2;
    x = [x1-a, x1+a, x2-a,x2+a];
    y = [y1-a, y1+a, y2-a,y2+a];
    inner_domains_coordinates = [x(1:2),y(1:2);...
        x(3:4),y(3:4)];
end

nu = 2;
nv = 2;
ConPts = zeros(n_nurbs_patch,nu,nv,DIM);
% The first parameter  is the number of NURBS patches (inner domains)
% The second parameter is the number of basis functions in x-direction
% The third parameter is the number of basis functions in y-direction
% The fourth parameter is the number spatial dimensions
for  e = 1:n_nurbs_patch
    x = inner_domains_coordinates(e,1:2);
    y = inner_domains_coordinates(e,3:4);
    ConPts(e,:,:,1) = [x(1) x(1)  ; x(2) x(2)];
    ConPts(e,:,:,2) = [y(1) y(2)  ; y(1) y(2)];
end
pu = 1; pv = 1; % The ultimate degree of B splines basis functions used is (pu+t).
weights =[1 1;1 1];
knotU   =[0 0  1 1];
knotV   =[0 0  1 1];

nurbs_original = cell(n_nurbs_patch,1);
for e = 1:n_nurbs_patch
    nurbs_original{e}.ConPts  = zeros(nu,nv,DIM);
    nurbs_original{e}.weights = weights;
    nurbs_original{e}.pu      = pu;
    nurbs_original{e}.pv      = pv;
    nurbs_original{e}.knotU   = knotU;
    nurbs_original{e}.knotV   = knotV;
    for i = 1:nu
        for j=1:nv
            for d = 1:DIM
                nurbs_original{e}.ConPts(i,j,d) = ConPts(e,i,j,d);
            end
        end
    end
end

nurbs_refine_domains = cell(n_nurbs_patch,1);
for e = 1:n_nurbs_patch
    [knotU,knotV]  = IGADegreeElevSurface(knotU,knotV,t); % Elevate both parametric degrees by the requested order.
    pu = nurbs_original{e}.pu + t;     pv = nurbs_original{e}.pv + t;
    nurbs_refine_domains{e} = IGA_2D_Grid(knotU,knotV,pu,pv, Refinement, Example);
end

%% Accumulate global degree-of-freedom offsets across patches
n_dofs = nurbs_refine_domains{1}.n_dofs_domains;
for e = 1:(n_nurbs_patch-1)
    nurbs_refine_domains{e+1}.n_dofs_domains  = nurbs_refine_domains{e+1}.n_dofs_domains + n_dofs;
    n_dofs =  nurbs_refine_domains{e+1}.n_dofs_domains;
end

for e = 2:n_nurbs_patch % Merge degrees of freedom shared by adjacent patches.
    nurbs_refine_domains{e} = Update_Edge_DoFs( nurbs_refine_domains{e-1}.n_dofs_domains, nurbs_refine_domains{e}   );
end

%%  Current uh
uh_filename   =  strcat('./Numerical Results/EigenFunction_Data/',Example,'/Nc_',num2str(Nc), '/p_',strcat(num2str(pu)));
uh_filename   =  strcat( strcat(uh_filename,'/uh_'),num2str(Refinement));
load(uh_filename, 'uh')
current_uh = uh;

%% Reference uh
if strcmp(Example,'Example_1')
    ref_nurbs_filename   =  strcat('./Numerical Results/EigenFunction_Data/',Example,'/Nc_25/','p_',strcat(num2str(3)));
    ref_nurbs_filename   =  strcat( strcat(ref_nurbs_filename,'/nurbs_'),num2str(7));
    ref_nurbs_filename   =  strcat(ref_nurbs_filename,'.mat');
    ref_uh_filename   =  strcat('./Numerical Results/EigenFunction_Data/',Example,'/Nc_25/','p_',strcat(num2str(3)));
    ref_uh_filename   =  strcat( strcat(ref_uh_filename,'/uh_'),num2str(7));
    ref_uh_filename   =  strcat(ref_uh_filename,'.mat');
elseif strcmp(Example,'Example_2')
    ref_nurbs_filename   =  strcat('./Numerical Results/EigenFunction_Data/',Example,'/Nc_30/','p_',strcat(num2str(3)));
    ref_nurbs_filename   =  strcat( strcat(ref_nurbs_filename,'/nurbs_'),num2str(5));
    ref_nurbs_filename   =  strcat(ref_nurbs_filename,'.mat');
    ref_uh_filename   =  strcat('./Numerical Results/EigenFunction_Data/',Example,'/Nc_30/','p_',strcat(num2str(3)));
    ref_uh_filename   =  strcat( strcat(ref_uh_filename,'/uh_'),num2str(5));
    ref_uh_filename   =  strcat(ref_uh_filename,'.mat');
elseif strcmp(Example,'Example_3')
    ref_nurbs_filename   =  strcat('./Numerical Results/EigenFunction_Data/',Example,'/Nc_25/','p_',strcat(num2str(3)));
    ref_nurbs_filename   =  strcat( strcat(ref_nurbs_filename,'/nurbs_'),num2str(7));
    ref_nurbs_filename   =  strcat(ref_nurbs_filename,'.mat');
    ref_uh_filename   =  strcat('./Numerical Results/EigenFunction_Data/',Example,'/Nc_25/','p_',strcat(num2str(3)));
    ref_uh_filename   =  strcat( strcat(ref_uh_filename,'/uh_'),num2str(7));
    ref_uh_filename   =  strcat(ref_uh_filename,'.mat');
end

load(ref_nurbs_filename, 'nurbs')
ref_nurbs = nurbs;
load(ref_uh_filename, 'uh')
ref_uh = uh;

L2_err_Max = zeros(n_nurbs_patch,n_eigenvalues);
L2_err_Min = zeros(n_nurbs_patch,n_eigenvalues);
L2_err     = zeros(n_nurbs_patch,n_eigenvalues);

if strcmp(Example, 'Example_1')
    uh_pw         =  current_uh(nurbs_refine_domains{1}.n_dofs+1:end,:);
    ref_uh_pw     =  ref_uh(ref_nurbs{1}.n_dofs+1:end,:);
end

if strcmp(Example, 'Example_2')
    uh_pw         =  current_uh(nurbs_refine_domains{2}.n_dofs_domains + 1:end,:);
    ref_uh_pw     =  ref_uh(ref_nurbs{2}.n_dofs_domains+1:end,:);
end

Omega_Area = L*L;
for i = 1:n_eigenvalues
    for e = 1:n_nurbs_patch
        if e == 1
            row    = 1:nurbs_refine_domains{1}.n_dofs;
            ref_row = 1:ref_nurbs{1}.n_dofs;
        else
            row = (1  + nurbs_refine_domains{e-1}.n_dofs_domains ): nurbs_refine_domains{e}.n_dofs_domains;
            ref_row = (1  + ref_nurbs{e-1}.n_dofs_domains ): ref_nurbs{e}.n_dofs_domains;
        end
        L2_err_Max(e,i) = Compute_L2_Error(nurbs_original{e},nurbs_refine_domains{e},current_uh(row,i),ref_nurbs{e},ref_uh(ref_row,i));
        L2_err_Min(e,i) = Compute_L2_Error(nurbs_original{e},nurbs_refine_domains{e},-current_uh(row,i),ref_nurbs{e},ref_uh(ref_row,i));
        L2_err(e,i)    = min(L2_err_Max(e,i), L2_err_Min(e,i));
    end

    dx = 1e-2; dy = dx;
    x = -L/2:dx:L/2;
    ele_area = dx*dy;
    y = x;
    [x,y] = meshgrid(x,y);
    nx = length(x);
    ny = length(x);
    L2_err_Max_pw = zeros(1,n_eigenvalues);
    L2_err_Min_pw = zeros(1,n_eigenvalues);
    L2_err_pw     = zeros(1,n_eigenvalues);
    for i1 = 1:nx-1
        for j1 = 1:ny-1
            if x(i1,j1) + dx/2 >=-a &&  x(i1,j1)+dx/2<=a &&  y(i1,j1)+dy/2>=-a &&  y(i1,j1)+dy/2<=a
                continue;
            end
            uh_pw_current_xy = uh_pw_at_xy(uh_pw(:,i),p,x(i1,j1)+dx/2,y(i1,j1)+dy/2,L)/sqrt(Omega_Area);
            ref_uh_pw_xy     = uh_pw_at_xy(ref_uh_pw(:,i),p,x(i1,j1)+dx/2,y(i1,j1)+dy/2,L)/sqrt(Omega_Area);
            L2_err_Max_pw(i) =  L2_err_Max_pw(i) + abs( uh_pw_current_xy - ref_uh_pw_xy ).^2*ele_area;
            L2_err_Min_pw(i) =  L2_err_Min_pw(i) + abs( uh_pw_current_xy - (-ref_uh_pw_xy) ).^2*ele_area;
        end
    end
    L2_err_pw(i) = min( L2_err_Max_pw(i),  L2_err_Min_pw(i)   );
end
L2_err_domains = zeros(1,n_eigenvalues);
for i=1:n_eigenvalues
    L2_err_domains(i) =  sqrt( sum( L2_err(:,i).^2 ) + L2_err_pw(i) );
end
L2_err = L2_err_domains ;
end
