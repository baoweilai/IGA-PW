function nurbsInfo = IGA_2D_Grid(knotU,knotV,pu,pv,Refinement)
%Build the 2-D IGA mesh data.

[Ubar,Vbar,dof]=IGAknotRefineSurface(knotU,knotV, pu,pv,Refinement);

nurbsInfo.Ubar=Ubar;
nurbsInfo.Vbar=Vbar;
nurbsInfo.dof=dof;

UBreaks =unique(Ubar);
VBreaks =unique(Vbar);

nurbsInfo.UBreaks=UBreaks;
nurbsInfo.VBreaks=VBreaks;

% ----- Uspan(end)=Uspan(end)-pu-1;Vspan(end)=Vspan(end)-pv-1;

m=length(Ubar)-pu-1;
n=length(Vbar)-pv-1;

nurbsInfo.m=m;
nurbsInfo.n=n;
nurbsInfo.n_dofs = m*n;

nurbsInfo.n_dofs_domains = m*n;


nurbsInfo.pu = pu;
nurbsInfo.pv = pv;


uNoEs = length(UBreaks)-1;
vNoEs=length(VBreaks)-1;
NoEs=uNoEs*vNoEs;

nurbsInfo.uNoEs=uNoEs;
nurbsInfo.vNoEs=vNoEs;
nurbsInfo.NoEs=NoEs;

Eledof=(pu+1)*(pv+1);
Element=zeros(NoEs,Eledof);
knotSpanIndex=zeros(NoEs,2);
Coordinate=zeros(NoEs,4);

Neighbour = zeros(uNoEs,4);


%%  For bottom edge

u_ele_dofs = pu+1;

bottom_edge_dofs     = zeros(uNoEs,2*u_ele_dofs); % The two layers dofs near the bottom boundary

bottom_edge_dofs_1st = zeros(uNoEs,u_ele_dofs);   % The first layer dofs lying on the bottom boundary

bottom_edge_node = zeros(uNoEs,2);


for e = 1:uNoEs
    i = findspan(Ubar,pu,UBreaks(e));
    bottom_edge_node(e,:) = [UBreaks(e),UBreaks(e+1)];

    k = 1;

    for j=1:2  % For the bottom boundary
        for i1=(i-pu):i
            bottom_edge_dofs(e,k) = i1 + (j-1)*m;
            k = k+1;
        end
    end

    j = 1;  % The index j corresponding to the bottom boundary
    k = 1;
    for i1 = (i-pu):i
        bottom_edge_dofs_1st(e,k) = i1+(j-1)*m;
        k = k+1;
    end


end


%%  For top edge

u_ele_dofs = pu+1;

top_edge_dofs = zeros(uNoEs,2*u_ele_dofs);   % [DOFs on edge, DOFs near edge]

top_edge_dofs_1st = zeros(uNoEs,u_ele_dofs); %  % The first layer dofs lying on the top boundary

top_edge_node = zeros(uNoEs,2);


for e = 1:uNoEs
    i = findspan(Ubar,pu,UBreaks(e));
    top_edge_node(e,:) = [UBreaks(e),UBreaks(e+1)];

    k = 1;


    for j = (n-1):n
        for i1 = (i-pu):i
            top_edge_dofs(e,k) = i1 + (j-1)*m;
            k = k+1;
        end
    end

    j = n; % The index j corresponding to the top boundary
    k = 1;
    for i1 = (i-pu):i
        top_edge_dofs_1st(e,k) = i1 + (j-1)*m;
        k = k + 1;
    end
end


%%  For left edge

v_ele_dofs = pv + 1;
left_edge_dofs     = zeros(vNoEs,2*v_ele_dofs); % [DOFs on edge, DOFs near edge]
left_edge_dofs_1st = zeros(vNoEs,v_ele_dofs);
left_edge_node = zeros(vNoEs,2);


for e = 1:vNoEs
    j = findspan(Vbar,pv,VBreaks(e));
    left_edge_node(e,:) = [VBreaks(e),VBreaks(e+1)];

    k = 1;


    for j1 = (j-pv):j
        for i=1:2
            left_edge_dofs(e,k) = i + (j1-1)*m;
            k = k +1;
        end
    end

    i = 1; % The index i corresponding to the left boundary
    k = 1;
    for j1 = (j-pv):j
        left_edge_dofs_1st(e,k) = i + (j1-1)*m;
        k = k+1;
    end


end


%%  For right edge

v_ele_dofs = pv + 1;
right_edge_dofs     = zeros(vNoEs,2*v_ele_dofs); % [DOFs on edge, DOFs near edge]
right_edge_dofs_1st = zeros(vNoEs,v_ele_dofs);
right_edge_node     = zeros(vNoEs,2);


for e = 1:vNoEs
    j = findspan(Vbar,pv,VBreaks(e));
    right_edge_node(e,:) = [VBreaks(e),VBreaks(e+1)];

    k = 1;

    for j1 = (j-pv):j
        for i=(m-1):m
            right_edge_dofs(e,k) = i + (j1-1)*m;
            k = k+1;
        end
    end

    i = m; % The index i corresponding to the right boundary
    k = 1;

    for j1 = (j-pv):j
        right_edge_dofs_1st(e,k) = i + (j1-1)*m;
        k = k+1;
    end

end


nurbsInfo.bottom_edge_dofs = bottom_edge_dofs;
nurbsInfo.right_edge_dofs  = right_edge_dofs;
nurbsInfo.top_edge_dofs    = top_edge_dofs;
nurbsInfo.left_edge_dofs   = left_edge_dofs;


nurbsInfo.bottom_edge_dofs_1st = bottom_edge_dofs_1st;
nurbsInfo.top_edge_dofs_1st    = top_edge_dofs_1st;
nurbsInfo.left_edge_dofs_1st   = left_edge_dofs_1st;
nurbsInfo.right_edge_dofs_1st  = right_edge_dofs_1st;


nurbsInfo.bottom_edge_node = bottom_edge_node;
nurbsInfo.right_edge_node  = right_edge_node;
nurbsInfo.top_edge_node    = top_edge_node;
nurbsInfo.left_edge_node   = left_edge_node;


for j1=1:vNoEs
    for i1=1:uNoEs
        row=zeros(1,Eledof);
        e=i1+(j1-1)*uNoEs;

        e_right =  e + 1;
        e_left  =  e - 1;
        e_down  =  e - uNoEs;
        e_up    =  e + uNoEs;

        Neighbour(e,:) = [e_left,e_right,e_down,e_up];

        Coordinate(e,:)=[UBreaks(i1:i1+1),VBreaks(j1:j1+1)];
        i=findspan(Ubar,pu,UBreaks(i1));
        j=findspan(Vbar,pv,VBreaks(j1));
        knotSpanIndex(e,:)=[i,j];
        for k=0:pv
            temp=(k*(pu+1)+1):(k+1)*(pu+1);
            tmp=m*(j-pv-1+k)+(i-pu:i);
            row(temp)=tmp;
            % N_{1,1}, N_{2,1}, ... , N_{m,1}; N_{1,2},N_{2,2},..., N_{m,2}; ...;
            % N_{1,n},N_{2,n},..., N_{m,n}.
        end
        Element(e,:)=row;
    end
end


invalid_ele_idx =  Neighbour<1 |  Neighbour>NoEs;

Neighbour(invalid_ele_idx) = -1;

nurbsInfo.Neighbour = Neighbour;


nurbsInfo.Element=Element;
nurbsInfo.Coordinate=Coordinate;
nurbsInfo.knotSpanIndex=knotSpanIndex;


bottom_row = zeros(uNoEs,pu+1);
bottom_column=zeros(uNoEs,Eledof);

top_row =bottom_row;
top_column=bottom_column;


bottom_node=zeros(uNoEs,2);
bottom_span=zeros(uNoEs,1);

for e=1:uNoEs
    bottom_node(e,:)=[UBreaks(e),UBreaks(e+1)];

    i=findspan(Ubar,pu,UBreaks(e));
    bottom_span(e)=i;
    bottom_row(e,:)=i-pu:i;

    j=pv+1;
    for k=0:pv
        temp=(k*(pu+1)+1):(k+1)*(pu+1);
        tmp=m*(j-pv-1+k)+(i-pu:i);
        bottom_column(e,temp)=tmp;
    end

    top_row(e,:)=m*(n-1) + (i-pu:i);

    j=n;
    for k=0:pv
        temp=(k*(pu+1)+1):(k+1)*(pu+1);
        tmp=m*(j-pv-1+k)+(i-pu:i);
        top_column(e,temp)=tmp;
    end
end

nurbsInfo.bottom_node=bottom_node;
nurbsInfo.bottom_span=bottom_span;

nurbsInfo.bottom_row=bottom_row;
nurbsInfo.bottom_column=bottom_column;

nurbsInfo.top_row=top_row;
nurbsInfo.top_column=top_column;


left = zeros(vNoEs,pv+1);
right=left;

for e=1:vNoEs
    j=knotSpanIndex(e,2);
    i=j-pv:j;
    left(e,:)= 1 + (i-1)*m;
    right(e,:)=m + (i-1)*m;
end


nurbsInfo.left=left;
nurbsInfo.right=right;


bottom_dofs_2_layers =zeros(1,2*m); % The 2 layers of DOFs  near  v=0;
top_dofs_2_layers   = zeros(1,2*m);    % The 2 layers of DOFs  near  v=1;
left_dofs_2_layers   =zeros(1,2*n);
right_dofs_2_layers =zeros(1,2*n);

k=1;
for j=1:2
    for i=1:m
        index = i + (j-1)*m;
        bottom_dofs_2_layers(k) = index;
        k = k+1;
    end
end

k=1;
for j=(n-1):n
    for i=1:m
        index = i + (j-1)*m;
        top_dofs_2_layers(k) = index;
        k = k+1;
    end
end


k=1;
for j=1:n
    for i=1:2
        index = i + (j-1)*m;
        left_dofs_2_layers(k) = index;
        k = k+1;
    end
end

% left_dofs_2_layers'


k=1;
for j=1:n
    for i=(m-1):m
        index = i + (j-1)*m;
        right_dofs_2_layers(k) = index;
        k = k+1;
    end
end

nurbsInfo.bottom_dofs_2_layers = bottom_dofs_2_layers;
nurbsInfo.top_dofs_2_layers = top_dofs_2_layers;
nurbsInfo.left_dofs_2_layers = left_dofs_2_layers;
nurbsInfo.right_dofs_2_layers  = right_dofs_2_layers;


% The DOF index for (i,j) is:   i + (j-1)*m

bottom_dofs =1:m; % this is for v=0;
top_dofs   = m*(n-1) +(1:m);
left_dofs   =(0:(n-1))*m  +1;
right_dofs =(0:(n-1))*m  + m;

nurbsInfo.bottom_dofs = bottom_dofs;
nurbsInfo.n_dofs_bottom = length(bottom_dofs);

nurbsInfo.top_dofs = top_dofs;
nurbsInfo.n_dofs_top = length(top_dofs);


nurbsInfo.left_dofs   = left_dofs;
nurbsInfo.n_dofs_left = length(left_dofs);

nurbsInfo.right_dofs  = right_dofs;
nurbsInfo.n_dofs_right = length(right_dofs);

bottom_dofs_2nd_layer =m + (1:m);
top_dofs_2nd_layer   = m*(n-2) +(1:m);
left_dofs_2nd_layer   =(0:(n-1))*m  + 2;
right_dofs_2nd_layer =(0:(n-1))*m  + m-1;


nurbsInfo.bottom_dofs_2nd_layer = bottom_dofs_2nd_layer;
nurbsInfo.top_dofs_2nd_layer = top_dofs_2nd_layer;

nurbsInfo.left_dofs_2nd_layer = left_dofs_2nd_layer;
nurbsInfo.right_dofs_2nd_layer  = right_dofs_2nd_layer;
